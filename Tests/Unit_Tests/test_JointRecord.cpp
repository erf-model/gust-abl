// The joint run's record (JointRecord.H): written by the I/O rank, read back on every rank bit for bit, and refused
// when a line is missing or malformed.
#include "JointRecord.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {

// A directory for one test, removed afterwards, on the I/O rank.
struct ScopedDir
{
    std::string name;
    explicit ScopedDir (std::string n) : name(std::move(n))
    {
        if (amrex::ParallelDescriptor::IOProcessor()) { amrex::UtilCreateCleanDirectory(name, false); }
        amrex::ParallelDescriptor::Barrier();
    }
    ~ScopedDir ()
    {
        amrex::ParallelDescriptor::Barrier();
        if (amrex::ParallelDescriptor::IOProcessor()) { amrex::FileSystem::RemoveAll(name); }
        amrex::ParallelDescriptor::Barrier();
    }
};

void write_text (const std::string& dir, const std::string& text)
{
    if (amrex::ParallelDescriptor::IOProcessor()) { std::ofstream(dir + "/Joint") << text; }
    amrex::ParallelDescriptor::Barrier();
}

// 1 on every rank when the predicate holds on every rank (an EXPECT on another rank would not fail the test).
bool on_all_ranks (bool ok)
{
    int v = ok ? 1 : 0;
    amrex::ParallelDescriptor::ReduceIntMin(v);
    return v == 1;
}

} // namespace

TEST(JointRecord, RoundTripIsExactOnEveryRank)
{
    // Each joint feedback (nest.feedback = defect, relax, relax_consistent, inject).
    for (const char* fb : {"defect", "relax", "relax_consistent", "inject"}) {
        const ScopedDir d(std::string("ut_joint_record_") + fb);
        gust::JointRecord w;
        w.feedback = fb;
        w.exchanges = 150;
        // Not short decimals: a record printed with fewer digits than a hex float would not restore the floors exactly.
        w.outer_cfl_floor = amrex::Real(1.0) / amrex::Real(3.0);
        w.nests = {{"crest", amrex::Real(0.1), "outer_chk_crest00271"},
                   {"lee", amrex::Real(0), "runs/a b/outer_chk_lee00012"}};   // 0: no floor; a path with a space
        // The stop test's window: values a short decimal would round, and a 0.
        w.window = {amrex::Real(1.0) / amrex::Real(7.0), amrex::Real(0), amrex::Real(2.0e-9) / amrex::Real(3.0)};
        // The factor, not a short decimal; 0 with defect.
        w.factor = w.feedback == "defect" ? amrex::Real(0) : amrex::Real(1.0) / amrex::Real(3.0);
        gust::write_joint_record(d.name, w);

        gust::JointRecord r;
        const bool found = gust::read_joint_record(d.name, r);
        bool same = found && r.feedback == w.feedback && r.exchanges == w.exchanges
                    && r.outer_cfl_floor == w.outer_cfl_floor && r.nests.size() == w.nests.size();
        for (std::size_t n = 0; same && n < w.nests.size(); ++n) {
            same = r.nests[n].name == w.nests[n].name && r.nests[n].cfl_floor == w.nests[n].cfl_floor
                   && r.nests[n].chk == w.nests[n].chk;
        }
        same = same && r.has_feedback && r.window == w.window && r.has_factor && r.factor == w.factor;
        EXPECT_TRUE(on_all_ranks(same));
    }
}

TEST(JointRecord, VersionOneReadsWithoutAWindow)
{
    // A record from before the window line (GUST_JOINT 1): it reads, with no window and no feedback field.
    const ScopedDir d("ut_joint_record_v1");
    write_text(d.name, "GUST_JOINT 1\nfeedback relax\nexchanges 7\nouter 0x1p-2\nnest a 0x0p+0 chk_a00010\n");
    gust::JointRecord r;
    r.window = {amrex::Real(1)};
    r.has_feedback = true;
    EXPECT_TRUE(on_all_ranks(gust::read_joint_record(d.name, r) && !r.has_feedback && r.window.empty()
                             && r.nests.size() == 1 && r.nests[0].chk == "chk_a00010"));
}

TEST(JointRecord, VersionTwoReadsWithoutAFactor)
{
    // A record from before the factor line (GUST_JOINT 2): it reads, with a window and no factor.
    const ScopedDir d("ut_joint_record_v2");
    write_text(d.name, "GUST_JOINT 2\nfeedback relax\nexchanges 7\nouter 0x1p-2\nwindow 1 0x1p-30\n"
                       "nest a 0x0p+0 chk_a00010\n");
    gust::JointRecord r;
    r.has_factor = true;
    r.factor = amrex::Real(1);
    EXPECT_TRUE(on_all_ranks(gust::read_joint_record(d.name, r) && r.has_feedback && !r.has_factor && r.factor == 0
                             && r.window.size() == 1 && r.nests.size() == 1));
}

TEST(JointRecord, RestoresOnlyWithTheRecordedFeedbackAndFactor)
{
    gust::JointRecord r;
    r.feedback = "relax";
    r.has_feedback = true;
    r.has_factor = true;
    r.factor = amrex::Real(0.5);
    EXPECT_EQ(gust::joint_restore_refusal(r, "relax", amrex::Real(0.5)), "");
    EXPECT_EQ(gust::joint_restore_refusal(r, "relax", amrex::Real(0.75)),
              "its joint record holds the feedback of nest.feedback_factor = 0.5, not 0.75");
    // Close factors are two factors, printed as two: with the precision's decimal digits, and with all that
    // round-trip where those print them alike (the adjacent value).
    const bool is_double = sizeof(amrex::Real) == 8;
    if (is_double) {
        EXPECT_EQ(gust::joint_restore_refusal(r, "relax", amrex::Real(0.5000001)),
                  "its joint record holds the feedback of nest.feedback_factor = 0.5, not 0.5000001");
    }
    EXPECT_EQ(gust::joint_restore_refusal(r, "relax", std::nextafter(amrex::Real(0.5), amrex::Real(1))),
              std::string("its joint record holds the feedback of nest.feedback_factor = 0.5, not ")
                  + (is_double ? "0.50000000000000011" : "0.50000006"));
    EXPECT_EQ(gust::joint_restore_refusal(r, "inject", amrex::Real(0.5)),
              "its joint record holds the feedback of nest.feedback = relax");
    // Version 2 has no factor: any factor restores.
    r.has_factor = false;
    r.factor = amrex::Real(0);
    EXPECT_EQ(gust::joint_restore_refusal(r, "relax", amrex::Real(0.75)), "");
    EXPECT_EQ(gust::joint_restore_refusal(r, "inject", amrex::Real(0.75)),
              "its joint record holds the feedback of nest.feedback = relax");
    // defect has factor 0 on both sides.
    r.feedback = "defect";
    r.has_factor = true;
    EXPECT_EQ(gust::joint_restore_refusal(r, "defect", amrex::Real(0)), "");
    // Version 1 restores nothing.
    r.has_feedback = false;
    EXPECT_EQ(gust::joint_restore_refusal(r, "defect", amrex::Real(0)),
              "its joint record (GUST_JOINT 1) holds no feedback");
}

TEST(JointRecord, NoFileIsNotAJointCheckpoint)
{
    const ScopedDir d("ut_joint_record_none");
    gust::JointRecord r;
    r.exchanges = 9;
    EXPECT_TRUE(on_all_ranks(!gust::read_joint_record(d.name, r)));
}

TEST(JointRecord, LineEndsAreNotPartOfThePath)
{
    // A record edited on another system: \r\n line ends and trailing spaces do not join the checkpoint's path.
    const ScopedDir d("ut_joint_record_crlf");
    write_text(d.name, "GUST_JOINT 1\r\nfeedback relax\r\nexchanges 7\r\nouter 0x1p-2\r\n"
                       "nest a 0x0p+0 runs/a b/chk_a00010 \t\r\nnest b 0x1p-1 chk_b00010\r\n");
    gust::JointRecord r;
    EXPECT_TRUE(on_all_ranks(gust::read_joint_record(d.name, r) && r.feedback == "relax" && r.exchanges == 7
                             && r.nests.size() == 2 && r.nests[0].chk == "runs/a b/chk_a00010"
                             && r.nests[1].chk == "chk_b00010" && r.nests[1].cfl_floor == 0.5));
}

TEST(JointRecord, MalformedRecordsAbort)
{
    const ScopedDir d("ut_joint_record_bad");
    const std::string head = "GUST_JOINT 1\nfeedback defect\nexchanges 4\nouter 0x1p-2\n";
    const std::string head2 = "GUST_JOINT 2\nfeedback defect\nexchanges 4\nouter 0x1p-2\n";
    const std::string head3 = "GUST_JOINT 3\nfeedback relax\nexchanges 4\nouter 0x1p-2\nwindow 0\n";
    const std::string bad[] = {
        "GUST_JOINT 4\nfeedback defect\nexchanges 4\nouter 0x1p-2\nwindow 0\nfactor 0x0p+0\nnest a 0x0p+0 chk\n",  // another version
        head3 + "nest a 0x0p+0 chk\n",                                                    // version 3, no factor line
        head3 + "factor -0x1p-1\nnest a 0x0p+0 chk\n",                                   // a negative factor
        head3 + "factor nan\nnest a 0x0p+0 chk\n",                                       // not finite
        head3 + "factor 0x0p+0\nnest a 0x0p+0 chk\n",                                    // relax with no factor
        head3 + "factors 0x1p-1\nnest a 0x0p+0 chk\n",                                   // not a factor line
        "GUST_JOINT 3\nfeedback defect\nexchanges 4\nouter 0x1p-2\nwindow 0\nfactor 0x1p-1\nnest a 0x0p+0 chk\n",  // defect has none
        "GUST_JOINT 0\nfeedback defect\nexchanges 4\nouter 0x1p-2\nnest a 0x0p+0 chk\n",   // no such version
        head2 + "nest a 0x0p+0 chk\n",                                                    // version 2, no window line
        head2 + "window 2 0x1p-30\nnest a 0x0p+0 chk\n",                                  // fewer values than counted
        head2 + "window 1 -0x1p-30\nnest a 0x0p+0 chk\n",                                 // a negative change
        head2 + "window 1 inf\nnest a 0x0p+0 chk\n",                                      // not finite
        head2 + "window -1\nnest a 0x0p+0 chk\n",                                         // a negative count
        head2 + "windows 0\nnest a 0x0p+0 chk\n",                                         // not a window line
        "GUST_JOINT 1\nfeedback none\nexchanges 4\nouter 0x1p-2\nnest a 0x0p+0 chk\n",     // not a joint feedback
        "GUST_JOINT 1\nfeedback defect\nexchanges 0\nouter 0x1p-2\nnest a 0x0p+0 chk\n",   // no exchange taken
        "GUST_JOINT 1\nfeedback defect\nexchanges 4\nouter -0x1p-2\nnest a 0x0p+0 chk\n",  // a negative floor
        "GUST_JOINT 1\nfeedback defect\nexchanges 4\nouter nan\nnest a 0x0p+0 chk\n",      // not finite
        head,                                                                              // no nest
        head + "nest a 0x0p+0\n",                                                          // no checkpoint
        head + "nest a 0x1p-3x chk\n",                                                     // trailing text in a floor
        head + "nests a 0x0p+0 chk\n",                                                     // not a nest line
    };
    for (const auto& text : bad) {
        write_text(d.name, text);
        gust::JointRecord r;
        bool threw = false;
        try { gust::read_joint_record(d.name, r); } catch (const std::runtime_error&) { threw = true; }
        EXPECT_TRUE(on_all_ranks(threw)) << text;
    }
    // The same head with one good nest line reads.
    write_text(d.name, head + "nest a 0x0p+0 chk\n");
    gust::JointRecord r;
    EXPECT_TRUE(on_all_ranks(gust::read_joint_record(d.name, r) && r.nests.size() == 1 && r.outer_cfl_floor == 0.25));
    // Version 2 with an empty window and with two changes reads.
    write_text(d.name, head2 + "window 0\nnest a 0x0p+0 chk\n");
    EXPECT_TRUE(on_all_ranks(gust::read_joint_record(d.name, r) && r.has_feedback && r.window.empty()));
    write_text(d.name, head2 + "window 2 0x1p-30 0x0p+0\nnest a 0x0p+0 chk\n");
    EXPECT_TRUE(on_all_ranks(gust::read_joint_record(d.name, r) && r.window.size() == 2
                             && r.window[0] == amrex::Real(0x1p-30) && r.window[1] == 0 && r.nests.size() == 1));
    // Version 3 with a factor reads.
    write_text(d.name, head3 + "factor 0x1p-1\nnest a 0x0p+0 chk\n");
    EXPECT_TRUE(on_all_ranks(gust::read_joint_record(d.name, r) && r.has_factor && r.factor == amrex::Real(0.5)));
}
