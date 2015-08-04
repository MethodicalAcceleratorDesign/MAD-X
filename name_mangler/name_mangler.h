#ifndef name_mangler_hh
#define name_mangler_hh

#ifdef __cplusplus
extern "C" {
#endif

#define NAME_MANGLER_BASE 4
#define NAME_MANGLER_SUFFIX 2

typedef void NameMangler_t;

extern NameMangler_t *NameMangler_init();
extern const char    *NameMangler_mangle(NameMangler_t *mangler, const char *str, char *dest );
extern const char    *NameMangler_demangle(NameMangler_t *mangler, const char *str, char *dest );
extern void           NameMangler_free(NameMangler_t *mangler );


#ifdef __cplusplus
};
#endif

#endif /* name_mangler_hh */
