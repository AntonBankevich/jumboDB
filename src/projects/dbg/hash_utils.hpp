//
// Created by anton on 7/24/20.
//

#pragma once
typedef unsigned __int128 htype128;

template<class Key>
struct alt_hasher {
    size_t operator()( const Key& k ) const;
};

template <>
struct alt_hasher<htype128>
{
    size_t operator()( const htype128& x ) const
    {
        return (size_t(x) * 31) ^ size_t(x >> 64u);
    }
};

inline std::ostream &operator<<(std::ostream &os, htype128 val) {
    std::vector<size_t> res;
    while(val != 0) {
        res.push_back(val%10);
        val /= 10;
    }
    for(auto it = res.rbegin(); it != res.rend(); ++it) {
        os << *it;
    }
    return os;
}

inline std::istream &operator>>(std::istream &is, htype128 &val) {
    val = 0;
    std::string tmp;
    is >> tmp;
    for(char c : tmp) {
        val = val * 10 + c - '0';
    }
    return is;
}

inline void writeHashs(std::ostream &os, const std::vector<htype128> &hash_list) {
    os << hash_list.size() << std::endl;
    for(htype128 h : hash_list) {
        size_t * tmp = reinterpret_cast<size_t*>(&h);
        os << *tmp << " " << *(tmp + 1) << std::endl;
    }
}

inline void readHashs(std::istream &is, std::vector<htype128> &hash_list) {
    size_t len;
    is >> len;
    size_t a[2];
    for(size_t i = 0; i < len; i++) {
        is >> a[0] >> a[1];
        htype128 * tmp = reinterpret_cast<htype128 *>(a);
        hash_list.push_back(*tmp);
    }
}

inline std::string decimal_string(htype128 n) {
    std::vector<char> res;
    while(n) {
        res.push_back('0' + n%10);
        n = n / 10;
    }
    return std::string(res.rbegin(), res.rend());
}

inline htype128 string128(const std::string &s) {
    htype128 res = 0;
    for(char c : s)
        res = res * 10 + c;
    return res;
}