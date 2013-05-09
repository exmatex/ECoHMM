#ifndef __YAML_OUTPUT_H
#define __YAML_OUTPUT_H

#include <stdio.h>

extern FILE* yamlFile;

void yamlBegin();
void yamlEnd();

void yamlAppInfo(FILE* file);

#endif
