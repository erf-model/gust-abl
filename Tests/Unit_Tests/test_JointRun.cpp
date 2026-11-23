// The joint run's log (JointRun.H): a domain name up to
// kJointLogNameWidth characters is printed as it is; a longer one keeps its first 12 and last 9 characters around
// "...", so the fixed-size lines keep their numbers and two names that differ in a suffix of up to 9 characters stay apart.
#include "JointRun.H"

#include <gtest/gtest.h>

#include <cstdio>
#include <string>

TEST(JointRun, LogNamesAreTruncatedToTheColumnWidth)
{
    ASSERT_EQ(gust::kJointLogNameWidth, 24);
    EXPECT_EQ(gust::joint_log_name("nest.hill"), "nest.hill");
    EXPECT_EQ(gust::joint_log_name("outer domain"), "outer domain");
    EXPECT_EQ(gust::joint_log_name(""), "");
    const std::string w24 = "nest.abcdefghijklmnopqrs";   // 24 characters: kept
    ASSERT_EQ(w24.size(), 24u);
    EXPECT_EQ(gust::joint_log_name(w24), w24);
    const std::string w25 = w24 + "t";                     // 25: truncated
    EXPECT_EQ(gust::joint_log_name(w25), "nest.abcdefg...lmnopqrst");
    const std::string a = "nest.ridge_upstream_far_field_01", b = "nest.ridge_upstream_far_field_02";
    EXPECT_EQ(gust::joint_log_name(a), "nest.ridge_u..._field_01");
    EXPECT_EQ(gust::joint_log_name(a).size(), 24u);
    EXPECT_NE(gust::joint_log_name(a), gust::joint_log_name(b));
    // The per-domain column of the run's last line (JointRun.cpp: a 160-character buffer) keeps its numbers with the
    // longest name a nest may have.
    const std::string longest(200, 'x');
    char col[160];
    const int n = std::snprintf(col, sizeof(col), "; %s %d iterations %.3f s", gust::joint_log_name(longest).c_str(),
                                2147483647, 1.0e9);
    EXPECT_LT(n, static_cast<int>(sizeof(col)));
}
