#include <algorithm>
#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <map>

#include "mad_6track_name_mangler.h"

class NameMangler {

  std::map<std::string, std::vector<std::string> > IDs; // the IDs, and their suffixes

  static std::string index2suffix(int index )
  {
    const int base = 'z' - 'a' + 1 + 'Z' - 'A' + 1;
    std::string suffix;
    for (int i=0; i < NAME_MANGLER_SUFFIX; i++) {
      int digit;
      if (index>0) {
	digit = index % base;
	index /= base;
      } else {
	digit = 0;
      }
      char digit_char = digit < (base>>1) ? ('a' + digit) : ('A' + digit - (base>>1));
      suffix = digit_char + suffix;
    }
    if (index>0) {
      std::cerr << "error: maximum number of IDs reached!\n";
    }
    return suffix;
  }

  static int suffix2index(const std::string &suffix )
  {
    if (suffix.length()>NAME_MANGLER_SUFFIX) {
      return -1;
    }
    const int base = 'z' - 'a' + 1 + 'Z' - 'A' + 1;
    int index = 0;
    for (int i=0; i < NAME_MANGLER_SUFFIX; i++) {
      int digit;
      if (suffix[i]>='a' && suffix[i]<='z') {
	digit = suffix[i] - 'a';
      } else if (suffix[i]>='A' && suffix[i]<='Z') {
	digit = suffix[i] - 'A' + (base>>1);
      } else {
	return -1;
      }
      index *= base;
      index += digit;
    }
    return index;
  }
  
public:
  
  std::string mangle(const std::string &str )
  {
    if (str.length()<NAME_MANGLER_BASE+NAME_MANGLER_SUFFIX)
      return str;

    std::string base = str.substr(0, NAME_MANGLER_BASE);
    std::string ending = str.substr(NAME_MANGLER_BASE);
    std::string suffix;
    
    std::vector<std::string> &endings = IDs[base];
    std::vector<std::string>::iterator ending_itr = std::find(endings.begin(), endings.end(), ending);
    if (ending_itr == endings.end()) {
      suffix = index2suffix(endings.size());
      endings.push_back(ending);
    } else {
      suffix = index2suffix(ending_itr - endings.begin());
    }
    return base + suffix;
  }

  std::string demangle(const std::string &str )
  {
    if (str.length()<NAME_MANGLER_BASE+NAME_MANGLER_SUFFIX)
      return str;

    std::string base = str.substr(0, NAME_MANGLER_BASE);
    std::string suffix = str.substr(NAME_MANGLER_BASE, NAME_MANGLER_SUFFIX);
    std::vector<std::string> &endings = IDs[base];
    int index = suffix2index(suffix);
    if (index == -1 || index >= int(endings.size())) {
      std::cerr << "error: incorrect element name '" << str << "'" << std::endl;
      return std::string();
    }
    return base + endings[index];
  }
  
};

NameMangler_t *NameMangler_init()
{
  return (void *)new class NameMangler;
}

void NameMangler_free(NameMangler_t *mangler )
{
  delete reinterpret_cast<class NameMangler*>(mangler);
}

const char *NameMangler_mangle(NameMangler_t *mangler, const char *str, char *dest )
{
  std::string mangled_str = reinterpret_cast<class NameMangler*>(mangler)->mangle(str);
  strcpy(dest, mangled_str.c_str());
  return dest;
}

const char *NameMangler_demangle(NameMangler_t *mangler, const char *str, char *dest )
{
  std::string demangled_str = reinterpret_cast<class NameMangler*>(mangler)->demangle(str);
  strcpy(dest, demangled_str.c_str());
  return dest;
}

/* 
** TEST & USAGE
*/

/*
int main()
{
  NameMangler_t *mangler;

  if ((mangler = NameMangler_init())) {

    char buf[255];

    std::cout << NameMangler_mangle(mangler, "Andrea", buf) << std::endl;
    std::cout << NameMangler_mangle(mangler, "Andreaa", buf) << std::endl;
    std::cout << NameMangler_mangle(mangler, "Andreaaa", buf) << std::endl;
    std::cout << NameMangler_mangle(mangler, "Andrea", buf) << std::endl;
    std::cout << NameMangler_mangle(mangler, "Andrea", buf) << std::endl;
    std::cout << NameMangler_mangle(mangler, "Andrea", buf) << std::endl;
    std::cout << NameMangler_mangle(mangler, "Andreaaa", buf) << std::endl;
    std::cout << NameMangler_mangle(mangler, "Andre2", buf) << std::endl;
    std::cout << NameMangler_mangle(mangler, "Ciao123", buf) << std::endl;
    std::cout << NameMangler_mangle(mangler, "Ciao1234", buf) << std::endl;
    std::cout << NameMangler_mangle(mangler, "Ciao12345", buf) << std::endl;
    std::cout << NameMangler_mangle(mangler, "Ciao1234", buf) << std::endl;
    std::cout << NameMangler_mangle(mangler, "Ciao12", buf) << std::endl;
    std::cout << NameMangler_mangle(mangler, "Ciao1", buf) << std::endl;
    std::cout << NameMangler_mangle(mangler, "Ciao", buf) << std::endl;

    for (int i=0; i<20*20; i++) {
      std::string new_str = std::string("PIFF") + char('a' + (i%20)) + char((i/20 + 'A'));
      NameMangler_mangle(mangler, new_str.c_str(), buf);
      NameMangler_demangle(mangler, buf, buf);
      if (new_str != std::string(buf)) {
	std::cerr << "A big big error occurred\n";
        break;
      }
    }
    
    NameMangler_free(mangler);
  }
  
  return 0;
}
*/
