#pragma once
/**
 * include/g2g/span.hpp
 * ====================
 * Portable `g2g::span` — a non-owning view over a contiguous range.
 *
 * G2G targets C++20, and `std::span` is the *only* C++20 standard-library
 * feature it relies on.  Unfortunately some HPC clusters still ship GCC 9
 * (e.g. Ubuntu 20.04's 9.4.0), whose libstdc++ does not provide `<span>` at
 * all — so a `#include <span>` aborts the build with
 *     fatal error: span: No such file or directory
 *
 * To stay buildable on those toolchains without giving up C++20 elsewhere,
 * every G2G header/source uses `g2g::span` instead of `std::span`:
 *
 *   • If the standard library provides `std::span` (GCC ≥ 10, Clang ≥ 7 with
 *     libc++, MSVC), `g2g::span` is a thin alias for it — zero overhead, full
 *     fidelity.
 *   • Otherwise, we provide a minimal dynamic-extent fallback implementing
 *     exactly the subset G2G uses:
 *       construction from (ptr, len) or a contiguous container; size();
 *       empty(); data(); operator[]; begin()/end(); front()/back(); subspan().
 *
 * Define G2G_SPAN_FORCE_FALLBACK to exercise the fallback on a modern stdlib
 * (used by the test build to verify the GCC-9 path on newer compilers).
 */

#if defined(__has_include)
#  if __has_include(<span>)
#    include <span>
#  endif
#endif

#if defined(__cpp_lib_span) && !defined(G2G_SPAN_FORCE_FALLBACK)

namespace g2g {
template <class T>
using span = std::span<T>;
}  // namespace g2g

#else  // ── fallback for pre-C++20 standard libraries (e.g. GCC 9.4) ──

#include <cstddef>
#include <iterator>     // std::data, std::size
#include <type_traits>

namespace g2g {

/// Minimal non-owning view over a contiguous sequence (dynamic extent only).
/// Mirrors the slice of std::span's interface that G2G actually uses.
template <class T>
class span {
public:
    using element_type   = T;
    using value_type     = std::remove_cv_t<T>;
    using size_type      = std::size_t;
    using pointer        = T*;
    using reference      = T&;
    using iterator       = T*;
    using const_iterator = const T*;

    constexpr span() noexcept : data_(nullptr), size_(0) {}
    constexpr span(pointer ptr, size_type count) noexcept
        : data_(ptr), size_(count) {}

    // Construct from a contiguous container (vector, array, ...) whose element
    // pointer is qualification-convertible to T*.  Lvalues bind here, as for
    // std::span's owning ranges.
    template <class Container,
              class = std::enable_if_t<
                  !std::is_same_v<std::remove_cv_t<Container>, span> &&
                  std::is_convertible_v<
                      std::remove_pointer_t<
                          decltype(std::data(std::declval<Container&>()))> (*)[],
                      T (*)[]>>>
    constexpr span(Container& c) noexcept
        : data_(std::data(c)),
          size_(static_cast<size_type>(std::size(c))) {}

    // Temporaries are accepted only for a const element type, mirroring
    // std::span's range constructor (borrowed_range || is_const<element_type>):
    // `hungarian_match(sv(a), sv(b))` passes rvalue vectors to span<const T>.
    template <class Container,
              class = std::enable_if_t<
                  std::is_const_v<T> &&
                  !std::is_same_v<std::remove_cv_t<Container>, span> &&
                  std::is_convertible_v<
                      std::remove_pointer_t<
                          decltype(std::data(std::declval<const Container&>()))> (*)[],
                      T (*)[]>>>
    constexpr span(const Container&& c) noexcept
        : data_(std::data(c)),
          size_(static_cast<size_type>(std::size(c))) {}

    constexpr pointer   data()  const noexcept { return data_; }
    constexpr size_type size()  const noexcept { return size_; }
    constexpr bool      empty() const noexcept { return size_ == 0; }

    constexpr reference operator[](size_type i) const { return data_[i]; }
    constexpr reference front() const { return data_[0]; }
    constexpr reference back()  const { return data_[size_ - 1]; }

    constexpr iterator begin() const noexcept { return data_; }
    constexpr iterator end()   const noexcept { return data_ + size_; }

    constexpr span subspan(size_type offset, size_type count) const {
        return span(data_ + offset, count);
    }
    constexpr span subspan(size_type offset) const {
        return span(data_ + offset, size_ - offset);
    }

private:
    pointer   data_;
    size_type size_;
};

}  // namespace g2g

#endif
