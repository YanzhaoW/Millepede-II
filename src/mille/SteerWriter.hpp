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

#pragma once

#include <cstdint>
#include <format>
#include <fstream>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

constexpr auto DEFAULT_MILLEPEDE_PARAMETER_FILE = "parameters.txt";
constexpr auto DEFAULT_STEER_FILENAME = "steer.txt";
constexpr auto DEFAULT_DATA_FILENAME = "data.bin";

namespace millepede
{
    enum class PedeMethod : uint8_t
    {
        inversion,
        diagonalization,
        fullGORES,
        sparseGMRES,
        cholesky,
        bandchooseby,
        HIP
    };
    class SteerWriter
    {
      public:
        SteerWriter() = default;
        void set_filepath(std::string_view filepath) { filepath_ = filepath; }
        void set_data_filepath(std::string_view filepath) { data_filepath_ = filepath; }
        void set_parameter_file(std::string_view filename) { parameter_file_ = filename; }
        void set_working_dir(std::string_view dir) { working_dir_ = dir; }

        void add_parameter_default(int par_id, const std::pair<float, float>& values);
        void add_method(PedeMethod method, const std::pair<float, float>& values);
        void add_other_options(std::vector<std::string> options);
        void write();

      private:
        std::map<PedeMethod, std::pair<float, float>> methods_;
        std::string filepath_ = DEFAULT_STEER_FILENAME;
        std::string parameter_file_ = DEFAULT_MILLEPEDE_PARAMETER_FILE;
        std::string data_filepath_ = DEFAULT_DATA_FILENAME;
        std::unordered_map<int, std::pair<float, float>> parameter_defaults_;
        std::vector<std::vector<std::string>> other_options_;
        std::string working_dir_ = ".";

        void write_parameter_defaults();
        void write_data_file(std::ofstream& ofile);
        void write_methods(std::ofstream& ofile);
        void write_others(std::ofstream& ofile);
    };

} // namespace millepede

template <>
class std::formatter<millepede::PedeMethod>
{
  public:
    static constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }
    template <typename FmtContent>
    constexpr auto format(const millepede::PedeMethod& method, FmtContent& ctn) const
    {
        auto str = [method]() -> std::string_view
        {
            using enum millepede::PedeMethod;
            switch (method)
            {
                case inversion:
                    return "inversion";
                case diagonalization:
                    return "diagonalization";
                case fullGORES:
                    return "fullGORES";
                case sparseGMRES:
                    return "sparseGMRES";
                case cholesky:
                    return "cholesky";
                case bandchooseby:
                    return "bandchooseby";
                case HIP:
                    return "HIP";
                default:
                    return "invalid";
            }
            return "invalid";
        }();
        return std::format_to(ctn.out(), "{}", str);
    }
};
