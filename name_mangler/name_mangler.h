#ifndef name_mangler_h
#define name_mangler_h

#ifdef __cplusplus
extern "C" {
#endif

#define NAME_MANGLER_BASE   14
#define NAME_MANGLER_SUFFIX 2

  typedef void NameMangler_t;

  extern NameMangler_t *NameMangler_init();
  extern void           NameMangler_free(NameMangler_t *mangler );
  extern const char    *NameMangler_mangle(NameMangler_t *mangler, const char *str, char *dest );
  extern const char    *NameMangler_demangle(NameMangler_t *mangler, const char *str, char *dest );

#ifdef __cplusplus
};
#endif

#endif /* name_mangler_h */
