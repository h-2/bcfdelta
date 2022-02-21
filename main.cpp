
#include <seqan3/argument_parser/all.hpp>

#include "decode.hpp"
#include "encode.hpp"

static_assert(bio::compression_traits<bio::compression_format::bgzf>::available);

int main(int argc, char ** argv)
{
    seqan3::argument_parser top_level_parser{
      "bcfdelta",
      argc,
      argv,
      seqan3::update_notifications::off,
      {"encode", "decode"}
    };

    top_level_parser.info.version           = version;
    top_level_parser.info.date              = date;
    top_level_parser.info.short_description = "Losslessly compress VCF and BCF files.";
    top_level_parser.info.synopsis.push_back("bcfdelta encode input_file[.vcf.gz|.bcf] output_file[.vcf.gz|.bcf]");
    top_level_parser.info.synopsis.push_back("bcfdelta decode input_file[.vcf.gz|.bcf] output_file[.vcf.gz|.bcf]");

    try
    {
        top_level_parser.parse();

        seqan3::argument_parser & sub_parser = top_level_parser.get_sub_parser();

        if (sub_parser.info.app_name == "bcfdelta-encode")
        {
            encode_options_t options = parse_encode_arguments(sub_parser);
            encode(options);
        }
        else if (sub_parser.info.app_name == "bcfdelta-decode")
        {
            decode_options_t options = parse_decode_arguments(sub_parser);
            decode(options);
        }
        else
        {
            std::cerr << "Unknown subcommand: " << sub_parser.info.app_name << '\n';
        }
    }
    catch (seqan3::argument_parser_error const & ext) // catch user errors
    {
        seqan3::debug_stream << "[bcfdelta] " << ext.what() << "\n"; // customise your error message
        return 1;
    }
}
