// SHA-256 against the FIPS 180-4 / NIST example vectors.
// The reference archives (test_Reference.cpp) rest on it.
#include "Sha256.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <string>

namespace {

TEST (Sha256, NistVectors)
{
    EXPECT_EQ(gust::sha256_hex(""),
              "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    EXPECT_EQ(gust::sha256_hex("abc"),
              "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    EXPECT_EQ(gust::sha256_hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
              "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
    EXPECT_EQ(gust::sha256_hex(std::string(1000000, 'a')),
              "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
}

// Padding boundaries: messages of 55, 56, 63, 64 and 65 bytes, fed whole and
// one byte at a time, give the same digest.
TEST (Sha256, ChunkingDoesNotMatter)
{
    for (const std::size_t n : {55u, 56u, 63u, 64u, 65u, 1000u}) {
        std::string msg(n, '\0');
        for (std::size_t i = 0; i < n; ++i) { msg[i] = static_cast<char>((i * 31 + 7) & 0xff); }
        gust::Sha256 s;
        for (const char c : msg) { s.update(&c, 1); }
        EXPECT_EQ(s.hex(), gust::sha256_hex(msg)) << n;
    }
}

// hex() leaves the object usable: asking twice gives the same digest, and
// more data afterwards gives the digest of the longer message.
TEST (Sha256, HexIsRepeatable)
{
    gust::Sha256 s;
    s.update("ab", 2);
    EXPECT_EQ(s.hex(), s.hex());
    s.update("c", 1);
    EXPECT_EQ(s.hex(), gust::sha256_hex("abc"));
}

TEST (Sha256, FileDigestMatchesBytes)
{
    const std::string name = "ut_sha256_" + std::to_string(amrex::ParallelDescriptor::MyProc()) + ".bin";
    std::string msg(200000, '\0');
    for (std::size_t i = 0; i < msg.size(); ++i) { msg[i] = static_cast<char>(i * 2654435761u >> 24); }
    {
        std::ofstream os(name, std::ios::binary);
        os.write(msg.data(), static_cast<std::streamsize>(msg.size()));
    }
    EXPECT_EQ(gust::sha256_file(name), gust::sha256_hex(msg));
    std::remove(name.c_str());
    EXPECT_EQ(gust::sha256_file(name), "");   // missing file
}

} // namespace

// The hex-float text that goes into the restart hashes (Buildings, Canopy) has the same bytes on every
// platform. MSVC pads "%a" to 13 hex digits; glibc and Apple's libc print the shortest form. The pins are
// the glibc forms, which the hashes pinned in test_BuildingForcing and test_BuildingTke were made with.
TEST (Sha256, HexfloatTextIsTheSameOnEveryPlatform)
{
    // MSVC's text, as it prints it, to the glibc form.
    EXPECT_EQ(gust::canonical_hexfloat_text("0x1.0000000000000p-1"), "0x1p-1");
    EXPECT_EQ(gust::canonical_hexfloat_text("0x1.99999999999a0p-4"), "0x1.99999999999ap-4");
    EXPECT_EQ(gust::canonical_hexfloat_text("-0x1.8000000000000p+1"), "-0x1.8p+1");
    EXPECT_EQ(gust::canonical_hexfloat_text("0x0.0000000000000p+0"), "0x0p+0");
    // Text already in that form, and text with no fraction, unchanged.
    EXPECT_EQ(gust::canonical_hexfloat_text("0x1.999999999999ap-4"), "0x1.999999999999ap-4");
    EXPECT_EQ(gust::canonical_hexfloat_text("0x1p-1"), "0x1p-1");
    EXPECT_EQ(gust::canonical_hexfloat_text("inf"), "inf");
    // This platform's printf, through the canonical form.
    EXPECT_EQ(gust::hexfloat_canonical(0.5), "0x1p-1");
    EXPECT_EQ(gust::hexfloat_canonical(0.1), "0x1.999999999999ap-4");
    EXPECT_EQ(gust::hexfloat_canonical(-3.0), "-0x1.8p+1");
    EXPECT_EQ(gust::hexfloat_canonical(0.0), "0x0p+0");
    EXPECT_EQ(gust::hexfloat_canonical(0.0001), "0x1.a36e2eb1c432dp-14");
}
