#pragma once

#include <concepts>
#include <functional>

#include <bio/var_io/header.hpp>

inline constexpr std::string_view version = "0.1.0";
inline constexpr std::string_view date    = "2022-02-18";

struct delta_error : std::runtime_error
{
    delta_error(auto const &... args) : std::runtime_error{(bio::detail::to_string(args) + ...)} {}
};

template <typename T1, typename T2>
concept compatible_alph = (std::same_as<T1, char> && std::same_as<T2, char>) ||
                          (!std::same_as<T1, char> && !std::same_as<T2, char> && std::integral<T1> &&
                           std::integral<T2>) ||
                          (std::same_as<T1, float> && std::same_as<T1, float>);

// function alias
auto & formulaG = bio::detail::vcf_gt_formula;

template <typename op_t = std::minus<>, bool skip_problematic = true>
struct delta_visitor
{
    std::string_view const id{};
    int32_t const          number{};
    size_t const           n_alts{};

    bio::var_io::header const * const hdr_ptr = nullptr;

    static constexpr auto op = bio::detail::overloaded{
      // floats are not substracted/added but XORed instead
      [](float & cur, float const last)
      {
          int32_t & cur_i = *(reinterpret_cast<int32_t *>(&cur));
          cur_i ^= *(reinterpret_cast<int32_t const *>(&last));
      },
      // integers
      []<typename cur_t, typename last_t>(cur_t & cur, last_t const last)
      {
          if constexpr (std::same_as<op_t, std::minus<>>)
              cur -=
                (cur == bio::var_io::missing_value<cur_t> || last == bio::var_io::missing_value<last_t>) ? 0 : last;
          else
              cur +=
                (cur == bio::var_io::missing_value<cur_t> || last == bio::var_io::missing_value<last_t>) ? 0 : last;
      }};

    template <typename last_rng_t, typename cur_rng_t>
    void operator()(last_rng_t & last_rng, cur_rng_t & cur_rng) const
    {
        using last_alph           = seqan3::range_innermost_value_t<last_rng_t>;
        constexpr size_t last_dim = seqan3::range_dimension_v<last_rng_t>;
        using cur_alph            = seqan3::range_innermost_value_t<cur_rng_t>;
        constexpr size_t cur_dim  = seqan3::range_dimension_v<cur_rng_t>;

        if constexpr (compatible_alph<last_alph, cur_alph> && last_dim == cur_dim)
        {
            size_t const n_sample = std::min(last_rng.size(), cur_rng.size());

            if (size_t n_sample_hdr = hdr_ptr->column_labels.size() - 9; n_sample > n_sample_hdr)
            {
                throw delta_error{"Current range has more entries (",
                                  n_sample,
                                  ") than there are samples in header (",
                                  n_sample_hdr,
                                  ")."};
            }

            constexpr auto error_or_not = [](auto &&... args)
            {
                if constexpr (!skip_problematic)
                    throw delta_error{std::forward<decltype(args)>(args)...};
            };

            if constexpr (cur_dim == 1)
            {
                if (number != 1)
                    throw delta_error{"wrong dimension"};

                for (size_t j = 0; j < n_sample; ++j)
                    op(cur_rng[j], last_rng[j]);
            }
            else if constexpr (cur_dim == 2)
            {
                if constexpr (!std::same_as<cur_alph, char>)
                {
                    // If the concats are the same length, directly subtract the concats!
                    // This does not perform checks on the inner lengths, so certain errors may not be detected
                    // NOTE that this is not a problem per se, because the reverse operation will do the same
                    if (std::span cur_concat = cur_rng.concat(), last_concat = last_rng.concat();
                        cur_concat.size() == last_concat.size() && cur_rng.size() == last_rng.size())
                    {
                        for (size_t i = 0; i < cur_concat.size(); ++i)
                            op(cur_concat[i], last_concat[i]);

                        return;
                    }
                }

                switch (number)
                {
                    case 0:
                        throw delta_error{"Genotype fields cannot be in FLAG state."};
                        break;
                    case 1:
                        throw delta_error{"wrong dimension"};
                        break;
                    case bio::var_io::header_number::dot:
                        if (n_alts == 1) // assuming that number is still same per record
                        {
                            for (size_t i = 0; i < n_sample; ++i)
                            {
                                if (last_rng[i].size() != cur_rng[i].size())
                                    continue; // since this is dot, we can't assume anything anyways

                                for (size_t j = 0; j < last_rng[i].size(); ++j)
                                    op(cur_rng[i][j], last_rng[i][j]);
                            }
                        }
                        else if (id == "PL3") // this is n_alts per one
                        {
                            for (size_t i = 0; i < n_sample; ++i)
                            {
                                if (last_rng[i].size() != 1)
                                    continue;

                                for (size_t j = 0; j < cur_rng[i].size(); ++j)
                                    op(cur_rng[i][j], last_rng[i][0]);
                            }
                        }
                        // else it cannot be compressed
                        break;
                    case bio::var_io::header_number::A:
                        for (size_t i = 0; i < n_sample; ++i)
                        {
                            if (last_rng[i].size() != 1)
                            {
                                error_or_not("Last range size: ", last_rng[i].size(), ". Expected: ", 1, ".");
                                continue;
                            }
                            if (cur_rng[i].size() != n_alts)
                            {
                                error_or_not("Current range size: ", cur_rng[i].size(), ". Expected: ", n_alts, ".");
                                continue;
                            }

                            for (size_t j = 0; j < n_alts; ++j)
                                op(cur_rng[i][j], last_rng[i][0]);
                        }
                        break;
                    case bio::var_io::header_number::R:
                        for (size_t i = 0; i < n_sample; ++i)
                        {
                            if (last_rng[i].size() != 2)
                            {
                                error_or_not("Last range size: ", last_rng[i].size(), ". Expected: ", 2, ".");
                                continue;
                            }
                            if (cur_rng[i].size() != n_alts + 1)
                            {
                                error_or_not("Current range size: ",
                                             cur_rng[i].size(),
                                             ". Expected: ",
                                             n_alts + 1,
                                             ".");
                                continue;
                            }

                            op(cur_rng[i][0], last_rng[i][0]);

                            for (size_t j = 1; j < n_alts + 1; ++j)
                                op(cur_rng[i][j], last_rng[i][1]);
                        }
                        break;
                    case bio::var_io::header_number::G:
                        {
                            size_t const inner_size = formulaG(n_alts, n_alts) + 1;

                            for (size_t i = 0; i < n_sample; ++i)
                            {
                                if (last_rng[i].size() != 3)
                                {
                                    error_or_not("Last range size: ", last_rng[i].size(), ". Expected: ", 3, ".");
                                    continue;
                                }
                                if (cur_rng[i].size() != inner_size)
                                {
                                    error_or_not("Current range size: ",
                                                 cur_rng[i].size(),
                                                 ". Expected: ",
                                                 inner_size,
                                                 ".");
                                    continue;
                                }

                                // [0, 0] mapped to first value
                                op(cur_rng[i][0], last_rng[i][0]);

                                // [0, k>=1] mapped to second
                                for (size_t k = 1; k <= n_alts; ++k)
                                    op(cur_rng[i][formulaG(0, k)], last_rng[i][1]);

                                // [j>=1, k>=1] mapped to third
                                for (size_t j = 1; j <= n_alts; ++j)
                                    for (size_t k = j; k <= n_alts; ++k)
                                        op(cur_rng[i][formulaG(j, k)], last_rng[i][2]);
                            }
                            break;
                        }
                    default: // any number > 1
                        assert(number > 1);
                        for (size_t i = 0; i < n_sample; ++i)
                        {
                            if (last_rng[i].size() != number)
                            {
                                error_or_not("Last range size: ", last_rng[i].size(), ". Expected: ", number, ".");
                                continue;
                            }
                            if (cur_rng[i].size() != number)
                            {
                                error_or_not("Current range size: ", cur_rng[i].size(), ". Expected: ", number, ".");
                                continue;
                            }

                            for (size_t j = 0; j < last_rng[i].size(); ++j)
                                op(cur_rng[i][j], last_rng[i][j]);
                        }
                        break;
                }
            }
            else // cur_dim == 3
            {
                throw delta_error{"Handling of vector-of-strings in Genotype not implemented"};
            }
        }
        else
        {
            throw std::runtime_error{std::string{"Incompatible types between records.\nFunction signature: "} +
                                     std::string{__PRETTY_FUNCTION__}};
        }
    }
};
