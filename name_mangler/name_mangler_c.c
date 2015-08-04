#include <stdio.h>

#include "name_mangler.h"

int main()
{
  NameMangler_t *mangler;
  if ((mangler = NameMangler_init())) {

    char buf[255];

    printf("%s\n", NameMangler_mangle(mangler, "Andrea", buf));
    printf("%s\n", NameMangler_mangle(mangler, "Andreaa", buf));
    printf("%s\n", NameMangler_mangle(mangler, "Andreaaa", buf));
    printf("%s\n", NameMangler_mangle(mangler, "Andrea", buf));
    printf("%s\n", NameMangler_mangle(mangler, "Andrea", buf));
    printf("%s\n", NameMangler_mangle(mangler, "Andrea", buf));
    printf("%s\n", NameMangler_mangle(mangler, "Andreaaa", buf));
    printf("%s\n", NameMangler_mangle(mangler, "Andre2", buf));
    printf("%s\n", NameMangler_mangle(mangler, "Ciao123", buf));
    printf("%s\n", NameMangler_mangle(mangler, "Ciao1234", buf));
    printf("%s\n", NameMangler_mangle(mangler, "Ciao12345", buf));
    printf("%s\n", NameMangler_mangle(mangler, "Ciao1234", buf));
    printf("%s\n", NameMangler_mangle(mangler, "Ciao12", buf));
    printf("%s\n", NameMangler_mangle(mangler, "Ciao1", buf));
    printf("%s\n", NameMangler_mangle(mangler, "Ciao", buf));

    char str[255];
    for (int i=0; i<20*20; i++) {
      sprintf(str, "PIFF%c%c", (char)('a' + (i%20)), (char)((i/20 + 'A')));
      printf("%i\t%s\n", i, NameMangler_mangle(mangler, str, buf));
    }

    while(fgets(buf, 255, stdin)) {
      printf("%s\n", NameMangler_demangle(mangler, buf, buf));
    }
  
    NameMangler_free(mangler);
  }
  
  return 0;
}

