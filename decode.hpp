#pragma once

#include <seqan3/argument_parser/all.hpp>

#include <bio/var_io/reader.hpp>
#include <bio/var_io/writer.hpp>

#include "shared.hpp"

struct decode_options_t
{
    std::filesystem::path input;
    std::filesystem::path output;
    size_t                threads = std::max<size_t>(1, std::min<size_t>(8, std::thread::hardware_concurrency()));
};

decode_options_t parse_decode_arguments(seqan3::argument_parser & parser)
{
    parser.info.short_description = "Losslessly compress VCF and BCF files (decompression sub-program).";
    parser.info.version           = version;
    parser.info.date              = date;
    parser.info.synopsis.push_back("bcfdelta decode input_file[.vcf.gz|.bcf] output_file[.vcf.gz|.bcf]");

    decode_options_t options{};

    parser.add_positional_option(options.input,
                                 "The input file.",
                                 seqan3::input_file_validator{
                                   {"vcf", "vcf.gz", "bcf"}
    });
    parser.add_positional_option(options.output,
                                 "The output file."); //, seqan3::output_file_validator{{"vcf", "vcf.gz", "bcf"}});

    parser.add_subsection("Performance:");

    parser.add_option(options.threads,
                      '@',
                      "threads",
                      "Maximum number of threads to use.",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{1u, std::thread::hardware_concurrency() * 2});

    parser.parse();

    return options;
}

void undo_delta(bio::var_io::default_record<> const & ref_record,
                bio::var_io::default_record<> &       record,
                bio::var_io::header const &           in_hdr,
                std::vector<int32_t> &                vec32_buffer,
                std::vector<std::vector<int32_t>> &   vecvec32_buffer)
{
    for (auto it = record.genotypes().begin(); it != record.genotypes().end(); ++it)
    {
        bio::var_io::header::format_t const & format = in_hdr.formats[in_hdr.string_to_format_pos().at(it->id)];

        if (!format.other_fields.contains("Encoding") || format.other_fields.at("Encoding") != "Delta")
            continue;

        for (auto lit = ref_record.genotypes().begin(); lit != ref_record.genotypes().end(); ++lit)
        {
            if (it->id == lit->id)
            {
                // current vector might be over int8 while values might become int32
                switch (bio::var_io::value_type_id{it->value.index()})
                {
                    case bio::var_io::value_type_id::int8:
                    case bio::var_io::value_type_id::int16:
                        {
                            auto fun = [&](auto & rng)
                            {
                                using rng_t = decltype(rng);
                                if constexpr (std::same_as<rng_t, std::vector<int8_t> &> ||
                                              std::same_as<rng_t, std::vector<int16_t> &>)
                                {
                                    bio::detail::sized_range_copy(rng, vec32_buffer);
                                }
                                else
                                {
                                    throw std::runtime_error{"Unreachable code reached."};
                                }
                            };

                            std::visit(fun, it->value);
                            it->value = std::move(vec32_buffer);
                            break;
                        }
                    case bio::var_io::value_type_id::vector_of_int8:
                    case bio::var_io::value_type_id::vector_of_int16:
                        {
                            auto fun = [&](auto & rng)
                            {
                                using rng_t = decltype(rng);
                                if constexpr (std::same_as<rng_t, std::vector<std::vector<int8_t>> &> ||
                                              std::same_as<rng_t, std::vector<std::vector<int16_t>> &>)
                                {
                                    vecvec32_buffer.resize(std::ranges::size(rng));
                                    for (size_t i = 0; i < std::ranges::size(rng); ++i)
                                        bio::detail::sized_range_copy(rng[i], vecvec32_buffer[i]);
                                }
                                else
                                {
                                    throw std::runtime_error{"Unreachable code reached."};
                                }
                            };

                            std::visit(fun, it->value);
                            it->value = std::move(vecvec32_buffer);
                            break;
                        }
                    default:
                        break;
                }

                if (!bio::detail::type_id_is_compatible(bio::var_io::value_type_id{it->value.index()},
                                                        bio::var_io::value_type_id{lit->value.index()}))
                {
                    throw std::runtime_error{"Incompatible types in variants"};
                }

                delta_visitor<std::plus<>> visitor{it->id, format.number, record.alt().size(), &in_hdr};

                std::visit(visitor, lit->value, it->value);
                break;
            }
        }
    }
}

void decode(decode_options_t const & options)
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

    bio::var_io::writer writer{options.output};

    bio::var_io::header const & in_hdr  = reader.header();
    bio::var_io::header         out_hdr = in_hdr;

    if (!out_hdr.string_to_info_pos().contains("DELTA_COMP") || !out_hdr.string_to_info_pos().contains("DELTA_REF"))
    {
        std::cerr << "The input file does not seem to be delta-compressed already. Exiting.\n";
        std::exit(1);
    }

    /** clean up the out-header **/
    std::erase_if(out_hdr.infos,
                  [](bio::var_io::header::info_t const & info)
                  { return info.id == "DELTA_COMP" || info.id == "DELTA_REF"; });

    for (bio::var_io::header::format_t & format : out_hdr.formats)
    {
        if (format.other_fields.contains("Encoding"))
            format.other_fields.erase("Encoding");
    }

    writer.set_header(out_hdr);

    /** decode **/
    bio::var_io::default_record<>     ref_record;
    std::vector<int32_t>              vec32_buffer;
    std::vector<std::vector<int32_t>> vecvec32_buffer;

    // TODO add check that first record is REF
    for (bio::var_io::default_record<> & record : reader)
    {
        bool needs_decompression = false;
        bool is_reference        = false;

        for (bio::var_io::info_element<bio::ownership::deep> const & info : record.info())
        {
            if (info.id == "DELTA_REF")
                is_reference = true;

            if (info.id == "DELTA_COMP")
            {
                needs_decompression = true;
                if (record.alt().size() == 1) // multi-allelic are never reference
                    is_reference = true;
                break;
            }
        }

        std::erase_if(record.info(),
                      [](auto const & info) { return info.id == "DELTA_REF" || info.id == "DELTA_COMP"; });

        if (needs_decompression)
            undo_delta(ref_record, record, in_hdr, vec32_buffer, vecvec32_buffer);

        writer.push_back(record);

        if (is_reference)
        {
            // backup the record to be able to refer to it next iteration
            ref_record = std::move(record);
        }
    }
}
