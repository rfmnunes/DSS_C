#ifndef REGISTRY_INCLUDED
#define REGISTRY_INCLUDED


#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#define atoll _atoi64
#define strdup _strdup
#endif

//class reg_key_class 
//{
//public:
//	char *name;
//	char *value;
//
//	int type;   /* 0-string; 1-integer; 2-floating point */
//	
//	union {
//		long long int lli;
//		double dbl;
//	} kval;
//
//	reg_key_class *next;
//
//	// constructor
//
//	reg_key_class(char *parameter, section )   // check if parameter exists, if not, create it
//	{
//		//check for existance
//
//	
//
//		if (!search_parameter_key(parameter))// is null, add a new one
//		{
//			next = section->klist;
//			section->klist = this;
//			name = strdup(parameter);
//			value = NULL;
//			type = 0;
//			kval.lli = 0;
//
//		
//		}
//		else 
//		{
//			value=NULL;
//		}
//
//
//
//		//if does not exist, add
//
//	};
//
//	//search for parameter on registy
//
//	reg_key_class* search_parameter_key(char *parameter)
//	{
//		if (strcmp(name, parameter))// 0 means equal
//		{
//			if (next)
//			{
//				return next->search_parameter_key(parameter);
//			}
//			else
//			{
//				return NULL;
//			}
//		}
//		else 
//		{
//			return this;
//		}
//
//		// else, return NULL
//
//
//	
//	};
//
//
//
//
//	//~reg_key_class()
//	//{
//	//	delete next;
//	//	if (name) free(name);
//	//	if (value) free(value);
//	//	delete k;
//	//};
//
//	///* parses the key as string and returns a pointer */
//	//char *get_string(reg_key *k) 
//	//{
//	//	k->type = 0;
//	//	return k->value;
//	//}; 
//	//
//	//	/* parses the key as an integer and returns the value */
//	//int get_int(reg_key *k) 
//	//{
//	//	return (int)get_llint(k);
//	//};
//
//	///* parses the key as a long integer and returns the value */
//	//long get_long(reg_key *k) 
//	//{
//	//	return (long)get_llint(k);
//	//};
//
//	///* parses the key as a long long integer and returns the value */
//	//long long int get_llint(reg_key *k) 
//	//{
//	//	if (k->type != 1) {
//	//		//k->kval.lli = strtoll(k->value, (char **) NULL, 10);      /* parse as integer */
//	//		k->kval.lli = atoll(k->value);
//	//		k->type = 1;
//	//	}
//	//	return k->kval.lli;
//	//};
//
//	///* parses the key as a float and returns the value */
//	//float get_float(reg_key *k) 
//	//{
//	//	return (float)get_double(k);
//	//};
//
//	///* parses the key as a double and returns the value */
//	//double get_double(reg_key *k) 
//	//{
//	//	if (k->type != 2) {
//	//		k->kval.dbl = atof(k->value);     /* parse as floating point */
//	//		k->type = 2;
//	//	}
//	//	return k->kval.dbl;
//	//};
//
//};
//TO DO FINISH THIS
//class registry_class 
//{
//public:
//	char *name;
//	reg_key_class *klist;
//	registry_class *next;
//
//
//	//constructor
//
//	registry_class(char *filename);//new registry
//	
//	registry_class(registry_class **r, char *filename);//merge constructor
//
//	// destructor
//	~registry_class() 
//	{
//		delete_registry(next);
//		delete_keylist(klist);
//		if (name) free(name);
//		free(this);
//	}
//	};
//
//	// methods
//	void dump_registry(registry_class *reg, char *filename) 
//	{
//		FILE *f;
//		reg_key *key;
//
//		f = fopen(filename, "w");
//		if (f) {
//			fprintf(f, "# Registry dump to file %s\n\n", filename);
//			while (reg) {
//				fprintf(f, "[%s]\n", reg->name);
//				key = reg->klist;
//				while (key) {
//					fprintf(f, "%s = %s\n", key->name, key->value);
//					key = key->next;
//				} /* while */
//				fprintf(f, "\n\n");
//				reg = reg->next;
//			} /* while */
//		} /* if */
//		fclose(f);
//	};
//
//	/* searches for any given entry and returns the key */
//
//	reg_key_class *get_key(registry_class *reg, char *section, char *parm)
//	{
//		
//		if (strcmp(reg->name, section))
//			return get_key(reg->next, section, parm);
//		else
//			return get_key_from_klist(reg->klist, section, parm);
//		
//		printf("Failed to get key %s:%s from registry!\n", section, parm);
//		return NULL;
//	};
//
//	reg_key_class *get_key_from_klist(reg_key_class *key, char* section, char *parm)
//	{
//		if (key) {
//			if (strcmp(key->name, parm))
//				return get_key_from_klist(key->next, section, parm);
//			else
//				return key;
//		}
//		printf("Failed to get key %s:%s from registry!\n", section, parm);
//		return NULL;
//	}; /* get_key */
//
//
//};


/* struct that defines a 64 bit resolution registry key */
typedef struct reg_key_type {
    char *name;
    char *value;
    int type;   /* 0-string; 1-integer; 2-floating point */
    union {
        long long int lli;
        double dbl;
    } kval;
    struct reg_key_type *next;
} reg_key;

/* struct for registry lists */
typedef struct registry_type {
    char *name;
    reg_key *klist;
    struct registry_type *next;
} registry;

/* builds a new registry */
registry *new_registry(char *filename);

/* deletes a registry from memory */
void delete_registry(registry *reg);

/* merges a file to given registry */
int merge_registry(registry **reg, char *filename);

/* searches for an given entry and returns the key */
reg_key *get_key(registry *reg, char *section, char *parm);

/* parses the key as string and returns a pointer */
char *get_string(reg_key *k);

/* parses the key as an integer and returns the value */
int get_int(reg_key *k);

/* parses the key as a long integer and returns the value */
long get_long(reg_key *k);

/* parses the key as a long long integer and returns the value */
long long int get_llint(reg_key *k);

/* parses the key as a float and returns the value */
float get_float(reg_key *k);

/* parses the key as a double and returns the value */
double get_double(reg_key *k);

void dump_registry(registry *r, char *filename);



int read_int_from_reg(registry *reg, char *group, char *parm);

#endif