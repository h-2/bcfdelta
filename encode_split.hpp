#pragma once

#include <bio/var_io/reader.hpp>

#include "shared.hpp"

struct split_buffers_t
{
    std::vector<int32_t>                                 ad_ref;
    seqan3::concatenated_sequences<std::vector<int32_t>> ad_alt;
    std::vector<int32_t>                                 pl1;
    seqan3::concatenated_sequences<std::vector<int32_t>> pl2;
    seqan3::concatenated_sequences<std::vector<int32_t>> pl3;
};

void do_split(bio::var_io::default_record<> & record, split_buffers_t & split_buffers)
{
    using genotype_t = bio::var_io::genotype_element<bio::ownership::deep>;

    size_t const n_alts  = record.alt().size();
    size_t const ad_size = n_alts + 1;
    size_t const pl_size = formulaG(n_alts, n_alts) + 1;

    for (auto it = record.genotypes().begin(); it != record.genotypes().end(); ++it)
    {
        if (it->id == "AD")
        {
            auto split_ad = bio::detail::overloaded{
              [&](auto & source) { throw delta_error{"AD field was not a collection of Integers."}; },
              [&]<typename int_t>(seqan3::concatenated_sequences<std::vector<int_t>> & source)
              {
                  genotype_t ad_ref{.id = "AD_REF", .value = std::move(split_buffers.ad_ref)};
                  auto &     ad_ref_vec = std::get<std::vector<int32_t>>(ad_ref.value);
                  ad_ref_vec.clear();
                  ad_ref_vec.reserve(source.size());

                  genotype_t ad_alt{.id = "AD_ALT", .value = std::move(split_buffers.ad_alt)};
                  auto &     ad_alt_vec = std::get<seqan3::concatenated_sequences<std::vector<int32_t>>>(ad_alt.value);
                  ad_alt_vec.clear();
                  ad_alt_vec.reserve(source.size());
                  ad_alt_vec.concat_reserve(source.size() * n_alts);

                  bool fail = false;
                  for (auto && inner_vec : source)
                  {
                      if (inner_vec.size() == 1) // only ref
                      {
                          ad_ref_vec.back() = inner_vec[0];
                          continue;
                      }
                      else if (inner_vec.size() != ad_size) // something is wrong
                      {
                          fail = true; // failure means we just reatin the original, unsplit field
                          break;
                      }

                      // first element goes to ad_ref
                      ad_ref_vec.push_back(inner_vec[0]);
                      // other elements go to ad_alt
                      ad_alt_vec.push_back(inner_vec.subspan(1));
                  }

                  if (!fail)
                  {
                      record.genotypes().erase(it); // this invalidates other iterators
                      record.genotypes().push_back(std::move(ad_ref));
                      record.genotypes().push_back(std::move(ad_alt));
                  }
              }};

            std::visit(split_ad, it->value);
            break;
        }
    }

    for (auto it = record.genotypes().begin(); it != record.genotypes().end(); ++it)
    {
        if (it->id == "PL")
        {
            auto split_pl = bio::detail::overloaded{
              [&](auto & source) { throw delta_error{"PL field was not a collection of Integers."}; },
              [&]<typename int_t>(seqan3::concatenated_sequences<std::vector<int_t>> & source)
              {
                  genotype_t pl1{.id = "PL1", .value = std::move(split_buffers.pl1)};
                  auto &     pl1_vec = std::get<std::vector<int32_t>>(pl1.value);
                  pl1_vec.reserve(source.size());

                  genotype_t pl2{.id = "PL2", .value = std::move(split_buffers.pl2)};
                  auto &     pl2_vec = std::get<seqan3::concatenated_sequences<std::vector<int32_t>>>(pl2.value);
                  pl2_vec.reserve(source.size());
                  pl2_vec.concat_reserve(source.size() * n_alts);

                  genotype_t pl3{.id = "PL3", .value = std::move(split_buffers.pl3)};
                  auto &     pl3_vec = std::get<seqan3::concatenated_sequences<std::vector<int32_t>>>(pl3.value);
                  //                            PL size - PL2-size - PL1-size
                  size_t     pl3_element_size = pl_size - n_alts - 1;
                  pl3_vec.reserve(source.size());
                  pl3_vec.concat_reserve(source.size() * pl3_element_size);

                  bool fail = false;

                  for (auto && inner_vec : source)
                  {
                      pl1_vec.emplace_back();
                      pl2_vec.push_back();
                      pl3_vec.push_back();

                      if (inner_vec.size() != pl_size) // something is wrong
                      {
                          if (inner_vec.size() == 0)
                          {
                              pl1_vec.back() = bio::var_io::missing_value<int32_t>;
                              continue; // empty vectors are OK/ignored
                          }
                          else
                          {
                              fail = true; // failure means we just reatin the original, unsplit field
                              break;
                          }
                      }

                      // [0, 0] mapped to first value
                      pl1_vec.back() = inner_vec[0];

                      // [0, k>=1] mapped to second
                      for (size_t k = 1; k <= n_alts; ++k)
                          pl2_vec.last_push_back(inner_vec[formulaG(0, k)]);

                      // [j>=1, k>=1] mapped to third
                      for (size_t j = 1; j <= n_alts; ++j)
                          for (size_t k = j; k <= n_alts; ++k)
                              pl3_vec.last_push_back(inner_vec[formulaG(j, k)]);
                  }

                  if (!fail)
                  {
                      record.genotypes().erase(it);
                      record.genotypes().push_back(std::move(pl1));
                      record.genotypes().push_back(std::move(pl2));
                      record.genotypes().push_back(std::move(pl3));
                  }
              }};

            std::visit(split_pl, it->value);
            break;
        }
    }
}

void salvage_split_buffers(bio::var_io::default_record<> & record, split_buffers_t & split_buffers)
{
    using vec_t     = std::vector<int32_t>;
    using vec_vec_t = seqan3::concatenated_sequences<std::vector<int32_t>>;

    for (auto it = record.genotypes().begin(); it != record.genotypes().end(); ++it)
    {
        if (it->id == "AD_REF")
        {
            auto & data = std::get<vec_t>(it->value);
            data.clear();
            split_buffers.ad_ref = std::move(data);
        }
        else if (it->id == "AD_ALT")
        {
            auto & data = std::get<vec_vec_t>(it->value);
            data.clear();
            split_buffers.ad_alt = std::move(data);
        }
        else if (it->id == "PL1")
        {
            auto & data = std::get<vec_t>(it->value);
            data.clear();
            split_buffers.pl1 = std::move(data);
        }
        else if (it->id == "PL2")
        {
            auto & data = std::get<vec_vec_t>(it->value);
            data.clear();
            split_buffers.pl2 = std::move(data);
        }
        else if (it->id == "PL3")
        {
            auto & data = std::get<vec_vec_t>(it->value);
            data.clear();
            split_buffers.pl3 = std::move(data);
        }
    }
}
