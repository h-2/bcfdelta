#pragma once

#include <seqan3/argument_parser/all.hpp>

#include <bio/var_io/reader.hpp>
#include <bio/var_io/writer.hpp>

#include "encode_delta.hpp"
#include "encode_split.hpp"
#include "shared.hpp"

struct encode_options_t
{
    std::filesystem::path input;
    std::filesystem::path output;
    uint64_t              ref_freq         = 10'000;
    bool                  delta_compress   = true;
    bool                  split_fields     = false;
    bool                  compress_ints    = true;
    bool                  compress_floats  = false;
    bool                  compress_chars   = false;
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
                                 seqan3::input_file_validator{
                                   {"vcf", "vcf.gz", "bcf"}
    });
    parser.add_positional_option(options.output,
                                 "The output file."); //, seqan3::output_file_validator{{"vcf", "vcf.gz", "bcf"}});

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

    parser.add_option(options.compress_ints, '\0', "compress-ints", "Delta-compress integers.");

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

void encode(encode_options_t const & options)
{
    size_t threads        = options.threads - 1; // subtract one for the main thread
    size_t reader_threads = threads / 3;
    size_t writer_threads = threads - reader_threads;

    auto reader_options =
      bio::var_io::reader_options{.field_types    = bio::var_io::field_types<bio::ownership::deep>,
                                  .stream_options = bio::transparent_istream_options{.threads = reader_threads + 1}};

    bio::var_io::reader reader{options.input, reader_options};

    auto writer_options =
      bio::var_io::writer_options{.stream_options = bio::transparent_ostream_options{.threads = writer_threads + 1}};

    bio::var_io::writer writer{options.output, writer_options};

    // "out_hdr" is a copy of "in_hdr"
    auto hdr = reader.header();

    if (options.split_fields)
    {
        // rename AD to AD_ALT
        hdr.formats.push_back({.id          = "AD_ALT",
                               .number      = bio::var_io::header_number::A,
                               .type        = "Integer",
                               .type_id     = bio::var_io::value_type_id::vector_of_int32,
                               .description = "ALT entries of AD field."});

        // add ad_ref
        hdr.formats.push_back({.id          = "AD_REF",
                               .number      = 1,
                               .type        = "Integer",
                               .type_id     = bio::var_io::value_type_id::int32,
                               .description = "REF entry of AD field."});

        hdr.formats.push_back({.id          = "PL1",
                               .number      = 1,
                               .type        = "Integer",
                               .type_id     = bio::var_io::value_type_id::int32,
                               .description = "PL values for 00."});

        hdr.formats.push_back({.id          = "PL2",
                               .number      = bio::var_io::header_number::A,
                               .type        = "Integer",
                               .type_id     = bio::var_io::value_type_id::vector_of_int32,
                               .description = "PL values for ab where a == 0 and b >= 1."});

        hdr.formats.push_back({.id          = "PL3",
                               .number      = bio::var_io::header_number::dot,
                               .type        = "Integer",
                               .type_id     = bio::var_io::value_type_id::vector_of_int32,
                               .description = "PL values for ab where a >= 1 and b >= 1"});
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
            bio::var_io::header::info_t info{.id          = "DELTA_COMP",
                                             .number      = 0,
                                             .type        = "Flag",
                                             .type_id     = bio::var_io::value_type_id::flag,
                                             .description = "Records with this flag have delta-compressed fields."};
            hdr.infos.push_back(std::move(info));
        }

        if (!hdr.string_to_info_pos().contains("DELTA_REF"))
        {
            bio::var_io::header::info_t info{.id      = "DELTA_REF",
                                             .number  = 0,
                                             .type    = "Flag",
                                             .type_id = bio::var_io::value_type_id::flag,
                                             .description =
                                               "This record is an 'anchor' for subsequent compressed records."};
            hdr.infos.push_back(std::move(info));
        }

        // all non-string fields are delta-compressed by default
        for (bio::var_io::header::format_t & format : hdr.formats)
        {
            bool do_compress = false;
            switch (format.type_id)
            {
                case bio::var_io::value_type_id::char8:
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

    split_buffers_t split_buffers;

    std::unique_ptr<bio::var_io::default_record<>> lrecord{new bio::var_io::default_record<>};
    std::unique_ptr<bio::var_io::default_record<>> brecord{new bio::var_io::default_record<>};
    lrecord->chrom() = "invalid";
    lrecord->pos()   = -1;

    for (bio::var_io::default_record<> & record : reader)
    {
        bio::var_io::default_record<> & last_record = *lrecord;
        bio::var_io::default_record<> & bak_record  = *brecord;

        /* split fields */
        if (options.split_fields)
            do_split(record, split_buffers);

        /* delta compression */
        if (options.delta_compress)
        {
            // backup the record as it is changed in-place
            bak_record = record;

            // this is a "reference record"
            if ((record.alt().size() == 1) && // multi-allelic can never be reference
                  (last_record.chrom() != record.chrom()) ||
                (last_record.pos() / options.ref_freq != record.pos() / options.ref_freq))
            {
                record.info().push_back({.id = "DELTA_REF", .value = true});
            }
            else // this will be delta-compressed
            {
                record.info().push_back({.id = "DELTA_COMP", .value = true});
                do_delta(last_record, record, hdr, options.skip_problematic);
            }
        }

        /* write the record */
        writer.push_back(record);

        /* get back some buffers */
        if (options.split_fields)
            salvage_split_buffers(record, split_buffers);

        /* make the backup of the current record the "last record" */
        if (options.delta_compress && record.alt().size() == 1)
            std::swap(lrecord, brecord);
    }
}
