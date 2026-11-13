#ifndef GUST_SHA256_H_
#define GUST_SHA256_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>

namespace gust {

//! SHA-256 (FIPS 180-4), for the checksums of reference archives. Host only.
class Sha256
{
public:
    Sha256 ();
    void update (const void* data, std::size_t n);
    //! The digest of everything given so far, as 64 lowercase hex
    //! characters. The object is unchanged and can take more data.
    [[nodiscard]] std::string hex () const;

private:
    std::string finish ();
    void block (const std::uint8_t* p);
    std::array<std::uint32_t, 8> m_h;
    std::array<std::uint8_t, 64> m_buf{};
    std::size_t m_buf_len{0};
    std::uint64_t m_total{0};
};

//! SHA-256 of a byte string, as hex.
std::string sha256_hex (const std::string& bytes);

//! SHA-256 of a file's contents, as hex. Returns an empty string if the file
//! cannot be opened or read.
std::string sha256_file (const std::string& path);

//! A value as hex-float text with the same bytes on every platform, for text that is hashed: printf's
//! "%a" with the mantissa's trailing zeros removed (MSVC pads it to 13 hex digits, glibc and Apple's libc
//! do not: 0.5 is "0x1.0000000000000p-1" on Windows and "0x1p-1" elsewhere).
std::string hexfloat_canonical (double v);

//! The canonical form of "%a" text: the trailing zeros of the mantissa's fraction removed, and its
//! point too when nothing is left after it. Text in that form already is returned unchanged.
std::string canonical_hexfloat_text (const std::string& a);

} // namespace gust

#endif
