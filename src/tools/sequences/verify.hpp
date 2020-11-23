//
// Created by anton on 19.12.2019.
//

#pragma once

#include <iostream>
#include <execinfo.h>
#include <cassert>

inline void print_stacktrace() {
    std::cout << "=== Stack Trace ===" << std::endl;

    const size_t max_stack_size = 1000;

    void *stack_pointers[max_stack_size];
    int count = backtrace(stack_pointers, max_stack_size);

    char **func_names = backtrace_symbols(stack_pointers, count);

    // Print the stack trace
    for (int i = 0; i < count; ++i)
        std::cerr << func_names[i] << std::endl;

    // Free the string pointers
    free(func_names);
}

#define VERIFY(expr)                                             \
    do {                                                         \
        if (!(expr)) {                                           \
            print_stacktrace();                                  \
            assert(expr);                                             \
            abort(); \
        };                                                       \
    } while(0);


