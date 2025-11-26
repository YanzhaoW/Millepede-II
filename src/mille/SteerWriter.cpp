/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2023 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "SteerWriter.hpp"
#include <filesystem>
#include <fstream>
#include <ios>
#include <ranges>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace millepede
{
    void SteerWriter::add_parameter_default(int par_id, const std::pair<float, float>& values)
    {
        if (auto iter = parameter_defaults_.find(par_id); iter != parameter_defaults_.end())
        {
            throw std::logic_error(std::format(
                "Default value of parameter {} has been already set to {}!", iter->second.first, iter->second.second));
        }
        parameter_defaults_.emplace(par_id, values);
    }

    void SteerWriter::add_method(PedeMethod method, const std::pair<float, float>& values)
    {
        if (auto iter = methods_.find(method); iter != methods_.end())
        {
            throw std::logic_error(std::format("Method {} has been already set to {}!", iter->first, iter->second));
        }
        methods_.emplace(method, values);
    }

    void SteerWriter::add_other_options(std::vector<std::string> options)
    {
        other_options_.push_back(std::move(options));
    }

    void SteerWriter::write()
    {
        auto ofile = std::ofstream{ std::filesystem::path{ working_dir_ } / filepath_,
                                    std::ios_base::out | std::ios_base::trunc };

        if (not ofile.is_open())
        {
            throw std::runtime_error(std::format("Failed to open file {}", filepath_));
        }

        write_parameter_defaults();
        write_data_file(ofile);
        write_methods(ofile);
        write_others(ofile);
        ofile << "end\n";
    }

    void SteerWriter::write_data_file(std::ofstream& ofile)
    {
        ofile << "Cfiles\n";
        ofile << data_filepath_ << "\n";
        ofile << parameter_file_ << "\n";
        ofile << "\n";
    }

    void SteerWriter::write_parameter_defaults()
    {
        auto ofile = std::ofstream{ std::filesystem::path{ working_dir_ } / parameter_file_,
                                    std::ios_base::out | std::ios_base::trunc };

        if (not ofile.is_open())
        {
            throw std::runtime_error(std::format("Can't open file {}", parameter_file_));
        }

        ofile << "Parameter\n";
        for (const auto& [par_id, values] : parameter_defaults_)
        {
            ofile << std::format("{} {:.3f} {:.3f}\n", par_id, values.first, values.second);
        }
    }

    void SteerWriter::write_methods(std::ofstream& ofile)
    {
        for (const auto& [method, values] : methods_)
        {
            ofile << std::format("method {} {} {}\n", method, values.first, values.second);
        }
        ofile << "\n";
    }

    void SteerWriter::write_others(std::ofstream& ofile)
    {
        for (const auto& other : other_options_)
        {
            ofile << std::format(
                "{}\n",
                std::views::join_with(other | std::views::transform([](const auto& val) -> std::string
                                                                    { return std::format("{}", val); }),
                                      " "));
        }
    }
} // namespace millepede
