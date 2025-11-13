#pragma once

#include <format>
#include <string>
#include <string_view>
#include <unordered_map>

namespace millepede
{
    struct ParResultEntry
    {
        int par_id = 0;
        float value = 0.F;
        float pre_sigma = 0.F;
        float value_diff = 0.F;
        float error = 0.F;
    };

    class ResultReader
    {
      public:
        ResultReader() = default;
        void set_filename(std::string_view filename) { filename_ = filename; }

        void read();
        void print();
        [[nodiscard]] auto get_pars() const -> const auto& { return par_results_; }

      private:
        std::string filename_;
        std::unordered_map<int, ParResultEntry> par_results_;
    };

} // namespace millepede

template <>
class std::formatter<millepede::ParResultEntry>
{
  public:
    static constexpr auto parse(format_parse_context& ctx) { return ctx.end(); }
    template <typename FmtContent>
    constexpr auto format(const millepede::ParResultEntry& entry, FmtContent& ctn) const
    {
        return std::format_to(ctn.out(),
                              "par id: {}, value: {}, sigma: {}, value_diff: {}, error: {}",
                              entry.par_id,
                              entry.value,
                              entry.pre_sigma,
                              entry.value_diff,
                              entry.error);
    }
};
