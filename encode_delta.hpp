#pragma once

#include <bio/var_io/reader.hpp>

#include "shared.hpp"

void do_delta(bio::var_io::default_record<> const & last_record,
              bio::var_io::default_record<> &       record,
              bio::var_io::header const &           hdr,
              bool const                            skip_problematic)
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
                    throw delta_error{"The type of this record's ",
                                      it->id,
                                      " field did is not compatible with the previous record."};
                }

                if (skip_problematic)
                {
                    delta_visitor<std::minus<>, true> visitor{it->id, format.number, record.alt().size(), &hdr};
                    std::visit(visitor, lit->value, it->value);
                }
                else
                {
                    delta_visitor<std::minus<>, false> visitor{it->id, format.number, record.alt().size(), &hdr};
                    std::visit(visitor, lit->value, it->value);
                }

                break;
            }
        }
    }
}
