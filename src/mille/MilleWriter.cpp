#include "MilleWriter.hpp"
#include "MilleEntry.hpp"
#include <cstddef>
#include <fstream>
#include <iostream>
#include <print>
#include <ranges>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace stdrng = std::ranges;

namespace
{
    template <typename T>
    inline auto get_size_in_bytes()
    {
        return static_cast<std::streamsize>(sizeof(T));
    }

    template <typename T>
    inline auto get_size_in_bytes(const std::vector<T>& data)
    {
        return static_cast<std::streamsize>(sizeof(T) * data.size());
    }

} // namespace

namespace millepede
{
    MilleWriter::MilleWriter(const std::string& out_filename, bool is_binary, bool is_write_zero)
        : is_binary_(is_binary)
        , is_zero_written_(is_write_zero)
    {
        output_file_.open(
            out_filename,
            (is_binary ? (std::ios::binary | std::ios::out | std::ios::trunc) : std::ios::out | std::ios::trunc));
        if (!output_file_.is_open())
        {
            throw std::runtime_error(std::format("Mille::Mille: Could not open {} as output file.", out_filename));
        }
    }

    void MilleWriter::mille(const MilleDataPoint& data_point)
    {
        if (data_point.sigma <= 0.)
        {
            return;
        }

        check_buffer_size(data_point.locals.size(), data_point.globals.size());

        if (buffer_.get_current_size() == 0)
        {
            buffer_.add_entry(0, 0.);
        }

        buffer_.add_entry(0, data_point.measurement);

        for (const auto [index, value] :
             stdrng::views::enumerate(data_point.locals) |
                 stdrng::views::filter([this](const auto& index_deriv) -> bool
                                       { return std::get<1>(index_deriv) != 0 or is_zero_written_; }))
        {
            buffer_.add_entry(static_cast<int>(index + 1), value);
        }

        buffer_.add_entry(0, data_point.sigma);

        for (const auto& [label, deriv] : data_point.globals)
        {
            if (deriv != 0 or is_zero_written_)
            {
                if ((label > 0 or is_zero_written_) and label <= max_label_size_)
                {
                    buffer_.add_entry(label, deriv);
                }
                else
                {
                    std::println(stderr, "Mille::mille: Invalid label {} <= 0 or > ", label);
                }
            }
        }
    }

    void MilleWriter::special(const std::vector<std::pair<int, float>>& special_data)
    {
        if (special_data.empty())
        {
            return;
        }
        if (has_special_done_)
        {
            throw std::logic_error("Mille::special: Special values already stored for this record.");
        }
        if (buffer_.get_current_size() == 0)
        {
            buffer_.add_entry(0, 0.);
        }

        buffer_.add_entry(0, 0.);
        buffer_.add_entry(0, -static_cast<float>(special_data.size()));
        for (const auto& [index, value] : special_data)
        {
            buffer_.add_entry(index, value);
        }
        has_special_done_ = true;
    }

    void MilleWriter::end()
    {
        if (buffer_.is_empty())
        {
            return;
        }
        is_binary_ ? write_to_binary() : write_to_non_binary();
        reset();
    }

    void MilleWriter::write_to_binary()
    {
        const auto data_size = static_cast<int>(buffer_.get_current_size());
        output_file_.write(reinterpret_cast<const char*>(&data_size), get_size_in_bytes<decltype(data_size)>());
        output_file_.write(reinterpret_cast<const char*>(buffer_.get_values().data()),
                           get_size_in_bytes(buffer_.get_values()));
        output_file_.write(reinterpret_cast<const char*>(buffer_.get_indices().data()),
                           get_size_in_bytes(buffer_.get_indices()));
    }

    void MilleWriter::write_to_non_binary()
    {
        output_file_ << buffer_.get_current_size() << "\n";
        output_file_ << std::format(
            "{}\n",
            std::views::join_with(buffer_.get_indices() | std::views::transform([](auto index) -> std::string
                                                                                { return std::to_string(index); }),
                                  ' '));
        output_file_ << std::format(
            "{}\n",
            std::views::join_with(buffer_.get_values() | std::views::transform([](auto value) -> std::string
                                                                               { return std::to_string(value); }),
                                  ' '));
    }

    void MilleWriter::reset()
    {
        buffer_.clear();
        has_special_done_ = false;
    }

    void MilleWriter::close() { output_file_.close(); }

    void MilleWriter::check_buffer_size(std::size_t n_local, std::size_t n_global)
    {
        if (buffer_.get_current_size() >= max_buffer_size_)
        {
            throw std::runtime_error(
                std::format("Mille::checkBufferSize: Buffer too short ({}), \n need space for nLocal ({}) \nGlobal "
                            "({}) local/global derivatives, {} already stored!",
                            max_buffer_size_,
                            n_local,
                            n_global,
                            buffer_.get_current_size()));
        }
    }
} // namespace millepede
