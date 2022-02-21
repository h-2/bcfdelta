#pragma once

#include <seqan3/argument_parser/all.hpp>

#include <bio/var_io/reader.hpp>
#include <bio/var_io/writer.hpp>

#include "shared.hpp"

struct encode_options_t
{
    std::filesystem::path input;
    std::filesystem::path output;
    uint64_t              ref_freq = 10'000;
    bool                  delta_compress = true;
    bool                  split_fields = false;
    bool                  compress_ints   = true;
    bool                  compress_floats = false;
    bool                  compress_chars  = false;
    bool                  skip_problematic = true;
    size_t                threads = std::max<size_t>(2, std::min<size_t>(8, std::thread::hardware_concurrency()));
};

encode_options_t parse_encode_arguments(seqan3::argument_parser & parser)
{
    parser.info.short_description = "Losslessly compress VCF and BCF files.";
    parser.info.version           = version;
    parser.info.date              = date;
    parser.info.synopsis.push_back("bcfdelta encode input_file[.vcf.gz|.bcf] output_file[.vcf.gz|.bcf]");

    encode_options_t options{};

    parser.add_positional_option(options.input,
                                 "The input file.",
                                 seqan3::input_file_validator{{"vcf", "vcf.gz", "bcf"}});
    parser.add_positional_option(options.output, "The output file.");//, seqan3::output_file_validator{{"vcf", "vcf.gz", "bcf"}});


    parser.add_subsection("Which data to compress:");

    parser.add_option(options.delta_compress,
                      'd',
                      "delta-compress",
                      "Encode genotype values as the difference to the previous record's values.",
                      seqan3::option_spec::hidden);

    parser.add_option(options.split_fields,
                      's',
                      "split-fields",
                      "Split certain fields so that their layout becomes better compressible.",
                      seqan3::option_spec::hidden);

    parser.add_option(options.compress_ints,
                      '\0',
                      "compress-ints",
                      "Delta-compress integers.");

    parser.add_option(options.compress_floats,
                      '\0',
                      "compress-floats",
                      "XOR-compress floats (Good for BCF output, possibly bad for VCF output).");

    parser.add_option(options.compress_chars,
                      '\0',
                      "compress-chars",
                      "Delta-compress characters (this does not refer to STRING fields, just to CHAR fields).");

    parser.add_option(options.skip_problematic,
                      '\0',
                      "skip-problematic",
                      "Skip sub-ranges that do not have expected size.");

    parser.add_subsection("Performance:");

    parser.add_option(options.threads,
                      '@',
                      "threads",
                      "Maximum number of threads to use.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{2u, std::thread::hardware_concurrency() * 2});

    parser.add_subsection("Tuning:");

    parser.add_option(options.ref_freq,
                      'f',
                      "ref-freq",
                      "Keep an uncompressed record every N basepairs.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{100, 1'000'000});

    parser.parse();

    return options;
}

void do_delta(bio::var_io::default_record<> const & last_record,
              bio::var_io::default_record<> &       record,
              bio::var_io::header const &           hdr,
              encode_options_t const &              options)
{
    for (auto it = record.genotypes().begin(); it != record.genotypes().end(); ++it)
    {
        bio::var_io::header::format_t const & format = hdr.formats[hdr.string_to_format_pos().at(it->id)];

        if (!format.other_fields.contains("Encoding") || format.other_fields.at("Encoding") != "Delta")
            continue;

        for (auto lit = last_record.genotypes().begin(); lit != last_record.genotypes().end(); ++lit)
        {
            if (it->id == lit->id)
            {
                if (!bio::detail::type_id_is_compatible(bio::var_io::value_type_id{it->value.index()},
                                                        bio::var_io::value_type_id{lit->value.index()}))
                {
                    throw delta_error{"The type of this record's ", it->id, " field did is not compatible with the previous record."};
                }

                if (options.skip_problematic)
                {
                    delta_visitor<std::minus<>, true> visitor{format.number, record.alt().size(), &hdr};
                    std::visit(visitor, lit->value, it->value);
                }
                else
                {
                    delta_visitor<std::minus<>, false> visitor{format.number, record.alt().size(), &hdr};
                    std::visit(visitor, lit->value, it->value);
                }

                break;
            }
        }
    }
}

void do_split(bio::var_io::default_record<> & record)
{
    for (auto it = record.genotypes().begin(); it != record.genotypes().end(); ++it)
    {
        if (it->id == "AD")
        {
            bio::var_io::genotype_element<bio::ownership::deep> ad_ref{.id = "AD_REF" };
            auto & ad_ref_vec = ad_ref.value.emplace<std::vector<int8_t>>();

            bio::var_io::genotype_element<bio::ownership::deep> ad_alt{.id = "AD_ALT" };;
            auto & ad_alt_vec = ad_alt.value.emplace<std::vector<std::vector<int8_t>>>();

            // the following assumes VCF input; TODO use std::visit
            auto & source = std::get<std::vector<std::vector<int32_t>>>(it->value);

            for (auto & inner_vec : source)
            {
                ad_ref_vec.push_back(inner_vec[0]);
                ad_alt_vec.emplace_back();
                std::ranges::copy(inner_vec.begin() + 1, inner_vec.end(), std::back_inserter(ad_alt_vec.back()));
            }

            record.genotypes().erase(it); // this invalidates other iterators
            record.genotypes().push_back(std::move(ad_ref));
            record.genotypes().push_back(std::move(ad_alt));
            break;
        }
    }

    for (auto it = record.genotypes().begin(); it != record.genotypes().end(); ++it)
    {
        if (it->id == "PL")
        {
            bio::var_io::genotype_element<bio::ownership::deep> pl_shared{.id = "PL_SHARED" };
            auto & pl_shared_vec = pl_shared.value.emplace<std::vector<std::vector<int16_t>>>();

            bio::var_io::genotype_element<bio::ownership::deep> pl_uniq{.id = "PL_UNIQ" };;
            auto & pl_uniq_vec = pl_uniq.value.emplace<std::vector<std::vector<int16_t>>>();

            // the following assumes VCF input; TODO use std::visit
            auto & source = std::get<std::vector<std::vector<int32_t>>>(it->value);

            for (auto & inner_vec : source)
            {
                pl_shared_vec.emplace_back();
                pl_uniq_vec.emplace_back();

                if (inner_vec.size() < 3)
                    break;

                //TODO properly deduce which values to put where
                pl_shared_vec.back().push_back(inner_vec[0]);
                pl_shared_vec.back().push_back(inner_vec[2]);
                pl_uniq_vec.back().push_back(inner_vec[1]);

                if (inner_vec.size() >= 6)
                {
                    pl_shared_vec.back().push_back(inner_vec[4]);
                    pl_shared_vec.back().push_back(inner_vec[5]);
                    pl_uniq_vec.back().push_back(inner_vec[3]);
                }

                if (inner_vec.size() > 6)
                    std::ranges::copy(inner_vec.begin() + 5, inner_vec.end(), std::back_inserter(pl_uniq_vec.back()));
            }

            record.genotypes().erase(it);
            record.genotypes().push_back(std::move(pl_shared));
            record.genotypes().push_back(std::move(pl_uniq));
            break;
        }
    }
}

void encode(encode_options_t const & options)
{
    size_t threads = options.threads  - 1; // subtract one for the main thread
    size_t reader_threads = threads / 3;
    size_t writer_threads = threads - reader_threads;

    auto reader_options = bio::var_io::reader_options{
        .field_types = bio::var_io::field_types<bio::ownership::deep>,
        .stream_options = bio::transparent_istream_options{ .threads = reader_threads + 1} };

    bio::var_io::reader reader{options.input, reader_options};

    auto writer_options = bio::var_io::writer_options{
        .stream_options = bio::transparent_ostream_options{ .threads = writer_threads + 1} };

    bio::var_io::writer writer{options.output, writer_options};

    // "out_hdr" is a copy of "in_hdr"
    auto hdr = reader.header();

    if (options.split_fields)
    {
        // rename AD to AD_ALT
        hdr.formats.push_back({.id          = "AD_ALT",
                               .number      = bio::var_io::header_number::A,
                               .type        = bio::var_io::value_type_id::vector_of_int8,
                               .description = "ALT entries of AD field."});
    //     hdr.formats.back().other_fields["Encoding"] = "RunLengthEncoding";

        // add ad_ref
        hdr.formats.push_back({ .id          = "AD_REF",
                                .number      = 1,
                                .type        = bio::var_io::value_type_id::int8,
                                .description = "REF entry of AD field."});

        hdr.formats.push_back({ .id          = "PL_SHARED",
                                .number      = bio::var_io::header_number::dot,
                                .type        = bio::var_io::value_type_id::int8,
                                .description = "PL values that are often shared in a record."});
    //     hdr.formats.back().other_fields["Encoding"] = "RunLengthEncoding";

        hdr.formats.push_back({ .id          = "PL_UNIQ",
                                .number      = bio::var_io::header_number::A,
                                .type        = bio::var_io::value_type_id::int8,
                                .description = "PL values that are usually specific to sample."});
    }

    if (options.delta_compress)
    {
        if (hdr.string_to_info_pos().contains("DELTA_COMP") || hdr.string_to_info_pos().contains("DELTA_REF"))
        {
            std::cerr << "The input file seems to be delta-compressed already. Exiting.\n";
            std::exit(1);
        }

        if (!hdr.string_to_info_pos().contains("DELTA_COMP"))
        {
            bio::var_io::header::info_t info{ .id           = "DELTA_COMP",
                                              .number       = 0,
                                              .type         = bio::var_io::value_type_id::flag,
                                              .description  = "Records with this flag have delta-compressed fields." };
            hdr.infos.push_back(std::move(info));
        }

        if (!hdr.string_to_info_pos().contains("DELTA_REF"))
        {
            bio::var_io::header::info_t info{ .id           = "DELTA_REF",
                                              .number       = 0,
                                              .type         = bio::var_io::value_type_id::flag,
                                              .description  = "This record is an 'anchor' for subsequent compressed records." };
            hdr.infos.push_back(std::move(info));
        }

        // all non-string fields are delta-compressed by default
        for (bio::var_io::header::format_t & format : hdr.formats)
        {
            bool do_compress = false;
            switch (format.type)
            {
                case bio::var_io::value_type_id::char8:
                case bio::var_io::value_type_id::vector_of_char8:
                    do_compress = options.compress_chars;
                    break;
                case bio::var_io::value_type_id::float32:
                case bio::var_io::value_type_id::vector_of_float32:
                    do_compress = options.compress_floats;
                    break;
                case bio::var_io::value_type_id::string:
                case bio::var_io::value_type_id::vector_of_string:
                    do_compress = false;
                    break;
                default: // integer cases
                    do_compress = options.compress_ints;
                    break;
            }

            if (do_compress)
                format.other_fields["Encoding"] = "Delta";
        }
    }

    writer.set_header(hdr);

    std::unique_ptr<bio::var_io::default_record<>> lrecord{new bio::var_io::default_record<>};
    std::unique_ptr<bio::var_io::default_record<>> brecord{new bio::var_io::default_record<>};
    lrecord->chrom() = "invalid";
    lrecord->pos() = -1;

    for (bio::var_io::default_record<> & record : reader)
    {
        bio::var_io::default_record<> & last_record = *lrecord;
        bio::var_io::default_record<> & bak_record  = *brecord;

        /* split fields */
        if (options.split_fields)
            do_split(record);

        /* delta compression */
        if (options.delta_compress)
        {
            // backup the record as it is changed in-place
            bak_record = record;

            // this is a "reference record"
            if ((record.alt().size() == 1) &&                   // multi-allelic can never be reference
                (last_record.chrom() != record.chrom()) ||
                (last_record.pos() / options.ref_freq != record.pos() / options.ref_freq))
            {
                record.info().push_back({.id = "DELTA_REF", .value = true});
            }
            else // this will be delta-compressed
            {
                record.info().push_back({.id = "DELTA_COMP", .value = true});
                do_delta(last_record, record, hdr, options);
            }
        }

        /* write the record */
        writer.push_back(record);

        if (options.delta_compress && record.alt().size() == 1)
            std::swap(lrecord, brecord);
    }
}
