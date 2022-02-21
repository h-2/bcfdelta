#pragma once

#include <concepts>
#include <functional>

#include <bio/var_io/header.hpp>

//!\brief Thrown if information given to output format didn't match expectations.
struct delta_error : std::runtime_error
{
    delta_error(auto const & ... args) : std::runtime_error{(bio::detail::to_string(args) + ...)}
    {}
};

template <typename T1, typename T2>
concept compatible_alph = (std::same_as<T1, char> && std::same_as<T2, char>) ||
                          (!std::same_as<T1, char> && !std::same_as<T2, char> && std::integral<T1> && std::integral<T2>) ||
                           (std::same_as<T1, float> && std::same_as<T1, float>);

template <typename op_t = std::minus<>, bool skip_problematic = true>
struct delta_visitor
{
    int32_t number{};
    size_t n_alts{};

    bio::var_io::header const * hdr_ptr = nullptr;

    static constexpr op_t op_impl{};

    static constexpr auto op = bio::detail::overloaded{
        [] (float const cur, float const last) -> float
        {
            union U
            {
                float f;
                int32_t i;
            };
            // undefined behaviour for the win!
            U u{cur};
            u.i ^= U{last}.i;
            return u.f;
        },
        [] <typename cur_t, typename last_t> (cur_t const cur, last_t const last) -> cur_t
        {
            return (cur == bio::var_io::missing_value<cur_t> || last == bio::var_io::missing_value<last_t>)
                    ? cur : op_impl(cur, last);
        }};

    template <typename last_rng_t, typename cur_rng_t>
    void operator()(last_rng_t & last_rng, cur_rng_t & cur_rng) const
    {
        using            last_alph = seqan3::range_innermost_value_t<last_rng_t>;
        constexpr size_t last_dim  = seqan3::range_dimension_v<last_rng_t>;
        using            cur_alph  = seqan3::range_innermost_value_t<cur_rng_t>;
        constexpr size_t cur_dim   = seqan3::range_dimension_v<cur_rng_t>;

        if constexpr (compatible_alph<last_alph, cur_alph> && last_dim == cur_dim)
        {
            size_t const n_sample = hdr_ptr->column_labels.size() - 9;
            if (last_rng.size() != n_sample)
                throw delta_error{"Last outer range size: ", last_rng.size(), ". Expected: ", n_sample, " (number of samples)."};
            if (cur_rng.size() != n_sample)
                throw delta_error{"Current outer range size: ", cur_rng.size(), ". Expected: ", n_sample, " (number of samples)."};

            constexpr auto error_or_not = [] (auto const & ... args)
            {
                if constexpr (!skip_problematic)
                    throw delta_error{args...};
            };

            if constexpr (cur_dim == 1)
            {
                if (number != 1)
                    throw delta_error{"wrong dimension"};

                for (size_t j = 0; j < last_rng.size(); ++j)
                    cur_rng[j] = op(cur_rng[j], last_rng[j]);
            }
            else if constexpr (cur_dim == 2)
            {
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
                            for (size_t i = 0; i < last_rng.size(); ++i)
                            {
                                if (last_rng[i].size() != cur_rng[i].size())
                                    continue; // since this is dot, we can't assume anything anyways

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
                                cur_rng[i][j] = op(cur_rng[i][j], last_rng[i][0]);
                        }
                        break;
                    case bio::var_io::header_number::R:
                        for (size_t i = 0; i < last_rng.size(); ++i)
                        {
                            if (last_rng[i].size() != 2)
                            {
                                error_or_not("Last range size: ", last_rng[i].size(), ". Expected: ", 2, ".");
                                continue;
                            }
                            if (cur_rng[i].size() != n_alts + 1)
                            {
                                error_or_not("Current range size: ", cur_rng[i].size(), ". Expected: ", n_alts + 1, ".");
                                continue;
                            }

                            cur_rng[i][0] = op(cur_rng[i][0], last_rng[i][0]);

                            for (size_t j = 1; j < n_alts + 1; ++j)
                                cur_rng[i][j] = op(cur_rng[i][j], last_rng[i][1]);
                        }
                        break;
                    case bio::var_io::header_number::G:
                    {
                        // see spec
                        auto formula = [] (size_t const a, size_t const b) { return (b*(b + 1)) / 2 + a; };

                        size_t const inner_size = formula(n_alts, n_alts) + 1;

                        for (size_t i = 0; i < last_rng.size(); ++i)
                        {
                            if (last_rng[i].size() != 3)
                            {
                                error_or_not("Last range size: ", last_rng[i].size(), ". Expected: ", 3, ".");
                                continue;
                            }
                            if (cur_rng[i].size() != inner_size)
                            {
                                error_or_not("Current range size: ", cur_rng[i].size(), ". Expected: ", inner_size, ".");
                                continue;
                            }

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
                        assert(number > 1);
                        for (size_t i = 0; i < last_rng.size(); ++i)
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
                                cur_rng[i][j] = op(cur_rng[i][j], last_rng[i][j]);
                        }
                        break;
                }
            }
            else
            {
                throw std::runtime_error{std::string{"Wrong dimensions between records.\nFunction signature:"} +
                                         std::string{__PRETTY_FUNCTION__}};
            }
        }
        else
        {
            throw std::runtime_error{std::string{"Incompatible types between records.\nFunction signature: "} +
                                     std::string{__PRETTY_FUNCTION__}};
        }
    }
};
