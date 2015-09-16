#ifndef NAME_MANGLER_H
#define NAME_MANGLER_H

#ifdef __cplusplus
extern "C" {
#endif

#define NAME_MANGLER_BASE   14
#define NAME_MANGLER_SUFFIX 2

const char *NameMangler_mangle(const char *str, char *dest );
const char *NameMangler_demangle(const char *str, char *dest );

#ifdef __cplusplus
}
#endif

#endif /* NAME_MANGLER_H */
