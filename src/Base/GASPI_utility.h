#ifndef SUCCESS_OR_DIE_H
#define SUCCESS_OR_DIE_H

#include <GASPI.h>
#include <stdlib.h>

void success_or_die (const char* file, const int line, const int ec);

#define SUCCESS_OR_DIE(ec) success_or_die(__FILE__, __LINE__, ec)
#endif
