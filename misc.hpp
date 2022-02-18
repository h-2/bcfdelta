#pragma once

#include <concepts>
#include <functional>

#include <bio/var_io/header.hpp>


template <typename T1, typename T2>
concept compatible_alph = (!std::same_as<T1, char> && !std::same_as<T2, char>) &&
                          ((std::integral<T1> && std::integral<T2>) ||
                           (std::same_as<T1, float> && std::same_as<T1, float>));

template <typename op_t = std::minus<>>
struct delta_visitor
{
    int32_t number{};
    size_t n_alts{};

    bio::var_io::header const * hdr_ptr = nullptr;

    op_t const op{};

    template <typename last_rng_t, typename cur_rng_t>
    void operator()(last_rng_t & last_rng, cur_rng_t & cur_rng) const
    {
        using            last_alph = seqan3::range_innermost_value_t<last_rng_t>;
        constexpr size_t last_dim  = seqan3::range_dimension_v<last_rng_t>;
        using            cur_alph  = seqan3::range_innermost_value_t<cur_rng_t>;
        constexpr size_t cur_dim   = seqan3::range_dimension_v<cur_rng_t>;

        if constexpr (compatible_alph<last_alph, cur_alph> && last_dim == cur_dim)
        {
            if (last_rng.size() != hdr_ptr->column_labels.size() - 9)
                throw std::runtime_error{"neq sample size"};
            if (cur_rng.size() != hdr_ptr->column_labels.size() - 9)
                throw std::runtime_error{"neq sample size"};

            if constexpr (cur_dim == 1)
            {
                if (number != 1)
                    throw std::runtime_error{"wrong dimension"};

                for (size_t j = 0; j < last_rng.size(); ++j)
                    cur_rng[j] = op(cur_rng[j], last_rng[j]);
            }
            else if constexpr (cur_dim == 2)
            {
                switch (number)
                {
                    case 0:
                        throw std::runtime_error{"Cannot have number 0 here"};
                        break;
                    case 1:
                        throw std::runtime_error{"wrong dimension"};
                        break;
                    case bio::var_io::header_number::dot:
                        if (n_alts == 1) // assuming that number is still same per record
                        {
                            for (size_t i = 0; i < last_rng.size(); ++i)
                            {
                                if (last_rng[i].size() != cur_rng[i].size())
                                    std::runtime_error{"inner length size"};

                                for (size_t j = 0; j < last_rng[i].size(); ++j)
                                    cur_rng[i][j] = op(cur_rng[i][j], last_rng[i][j]);
                            }

                        }
                        // else it cannot be compressed
                        break;
                    case bio::var_io::header_number::A:
                        for (size_t i = 0; i < last_rng.size(); ++i)
                        {
                            if (last_rng[i].size() != 1)
                                std::runtime_error{"miscount"};
                            if (cur_rng[i].size() != n_alts)
                                std::runtime_error{"miscount"};

                            for (size_t j = 0; j < n_alts; ++j)
                                cur_rng[i][j] = op(cur_rng[i][j], last_rng[i][0]);
                        }
                        break;
                    case bio::var_io::header_number::R:
                        for (size_t i = 0; i < last_rng.size(); ++i)
                        {
                            if (last_rng[i].size() != 2)
                                std::runtime_error{"miscount"};
                            if (cur_rng[i].size() != n_alts + 1)
                                std::runtime_error{"miscount"};

                            cur_rng[i][0] = op(cur_rng[i][0], last_rng[i][0]);

                            for (size_t j = 1; j < n_alts + 1; ++j)
                                cur_rng[i][j] = op(cur_rng[i][j], last_rng[i][1]);
                        }
                        break;
                    case bio::var_io::header_number::G:
                    {
                        // see spec
                        auto formula = [] (size_t const a, size_t const b) { return (b*(b + 1)) / 2 + a; };

                        size_t inner_size = formula(n_alts, n_alts) + 1;

                        for (size_t i = 0; i < last_rng.size(); ++i)
                        {
                            if (last_rng[i].size() != 3)
                                std::runtime_error{"miscount"};
                            if (cur_rng[i].size() != inner_size)
                                std::runtime_error{"miscount"};

                            // [0, 0] mapped to first value
                            cur_rng[i][0] = op(cur_rng[i][0], last_rng[i][0]);

                            // [0, k>=1] mapped to second
                            for (size_t k = 1; k <= n_alts; ++k)
                                cur_rng[i][formula(0, k)] = op(cur_rng[i][formula(0, k)], last_rng[i][1]);

                            // [j>=1, k>=1] mapped to third
                            for (size_t j = 1; j <= n_alts; ++j)
                                for (size_t k = j; k <= n_alts; ++k)
                                    cur_rng[i][formula(j, k)] = op(cur_rng[i][formula(j, k)], last_rng[i][2]);
                        }
                        break;
                    }
                    default: // any number > 1
                        for (size_t i = 0; i < last_rng.size(); ++i)
                        {
                            if (last_rng[i].size() != number || cur_rng[i].size() != number)
                                std::runtime_error{"inner length size"};

                            for (size_t j = 0; j < last_rng[i].size(); ++j)
                                cur_rng[i][j] = op(cur_rng[i][j], last_rng[i][j]);
                        }
                        break;
                }
            }
            else
            {
                throw std::runtime_error{"wrong dimension"};
            }
        }
        else
        {
            throw std::runtime_error{"Incompatible types between records."};
        }
    }
};
