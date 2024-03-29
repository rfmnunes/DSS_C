#include "registry.h"
#include <fstream>

// TO DO GET SPACES ON PATHS

/* macros that identify the registry parser states */
#define IDLE           0
#define COMMENT        1
#define SECTION        2
#define BEGIN_PARM     3
#define PARM           4
#define BEGIN_VALUE    5
#define VALUE          6
#define END_VALUE      7

/* registry prototypes */
int is_alpha(char c);
void delete_keylist(reg_key *k);
reg_key *get_key_from_klist(reg_key *key, char *section, char *parm);



/* starts a new registry */
registry *new_registry (char *filename) {
	registry *reg;
	int l;

	reg = NULL;
	if ((l = merge_registry(&reg, filename)) > 0) {
		delete_registry(reg);
		reg = NULL;
		printf("Registry parsing failed on file %s at line %d!\n", filename, l);
	}
	return reg;
} /* new_registry */

/* parses and merges a file to given registry */
int merge_registry(registry **r, char *filename)
{
	//FILE *f;
	std::ifstream f;
	char state, buf[1024];
	int c;
	int i, nl;
	registry *section, *reg;
	reg_key *key;

	reg = *r;
	section = reg = NULL;
	state = IDLE;
	i = 0;
	nl = 1;
	f.open(filename);

	//if ((f = fopen(filename, "r")) != NULL) {
	if (!f.fail()){
		c=f.get();
		//c = fgetc(f);
		while(c != EOF) {
			switch(state) {
			case IDLE:  /* parser waiting for something... */
				if (c == '#')
					state = COMMENT;
				else if (c == '[') {
					state = SECTION;
					i = 0;
				} else if (c == 10)  /* new line */
					nl++;
				else if (section && is_alpha(c)) {
					/* start parsing new parameter */
					i = 0;
					buf[i++] = c;
					state = PARM;
				}
				break;
			case COMMENT: /* parsing a comment */
				if (c == 10) {
					state = IDLE;
					nl++;
				}
				break;
			case SECTION: /* parsing a new section */
				if (c == ']') {
					/* finished parsing section ID */
					buf[i] = 0;
					//printf_dbg2("merge_registry(): section: >%s<\n", buf);
					/* search for section */
					section = reg;
					while (section)
						if (strcmp(section->name, buf))
							section = section->next;
						else
							break;
					/* if not found, add new section */
					if (!section) {
						section = (registry *) malloc(sizeof(registry));
						section->next = reg;
						reg = section;
						section->name = strdup(buf);
						section->klist = NULL;
					}
					state = BEGIN_PARM;
					i = 0;
				} else if (is_alpha(c))
					/* save char from section ID */
					buf[i++] = c;
				else
					return nl;
				break;
			case BEGIN_PARM:
				if (c == '#') {
					i = 0;
					state = COMMENT;
				} else if (is_alpha(c)) {
					buf[i++] = c;
					state = PARM;
				} else if (c == 10) {
					nl++;
				}
				break;
			case PARM: /* parsing a new parameter in a given section */
				if (c == '=') {
					buf[i] = 0;
					/* search parm */
					key = section->klist;
					while (key)
						if (strcmp(key->name, buf))
							key = key->next;
						else
							break;
					/* if not found, add new key to section */
					if (!key) {
						key = (reg_key *) malloc(sizeof(reg_key));
						key->next = section->klist;
						section->klist = key;
						key->name = strdup(buf);
						key->value = NULL;
						key->type = 0;
						key->kval.lli = 0;
					} else {
						free(key->value);
					}
					state = BEGIN_VALUE;
					i = 0;
				} else if (is_alpha(c))
					buf[i++] = c;
				else if ((c != ' ') && (c != '\t'))
					return nl;
				break;
			case BEGIN_VALUE:
				if (is_alpha(c)) {
					buf[i++] = c;
					state = VALUE;
				} else if ((c != ' ') && (c != '\t'))
					return nl;
				break;
			case VALUE: /* parsing a value for a given parameter */
				if (is_alpha(c)|| c == ' ' )
					buf[i++] = c;
				else {
					/* end of value */
					buf[i] = 0;
					/* store value */
					key->value = strdup(buf);
					//printf_dbg2("merge_registry(): parm: >%s< >%s<\n", key->name, key->value);
					state = END_VALUE;
					i = 0;
					/* value of c is going to get lost in loop, we have to check if its a newline here */
					if( c == 10) {
						state = IDLE;
						nl++;
					}
				}
				break;
			case END_VALUE:
				if (c == 10) {
					state = IDLE;
					nl++;
				}
				break;
			default:
				/* do nothing */
				break;
			} /* switch */
			/* get next char */
			//c = fgetc(f);
			c=f.get();
		} /* while */
		/* parse file */
	}
	//close file
	//fclose (f);
	f.close();
	*r = reg;
	return 0;
} /* merge_registry */

/* deletes a set of keys from memory */
void delete_keylist(reg_key *k) {
	if (k) {
		delete_keylist(k->next);
		if (k->name) free(k->name);
		if (k->value) free(k->value);
		free(k);
	}
} /* delete_keylist */



/* deletes the registry from memory */
void delete_registry(registry *r) {
	if (r) {
		delete_registry(r->next);
		delete_keylist(r->klist);
		if (r->name) free(r->name);
		free(r);
	}
} /* delete_registry */

/* searches for the key on a given section and returns the key */
reg_key *get_key_from_klist(reg_key *key, char* section, char *parm)
{
	if (key) {
		if (strcmp(key->name, parm))
			return get_key_from_klist(key->next, section, parm);
		else
			return key;
	}
	printf("Failed to get key %s:%s from registry!\n", section, parm);
	return NULL;
} /* get_key */

/* searches for any given entry and returns the key */
reg_key *get_key(registry *reg, char *section, char *parm)
{
	if (reg) {
		if (strcmp(reg->name, section))
			return get_key(reg->next, section, parm);
		else
			return get_key_from_klist(reg->klist, section, parm);
	}
	printf("Failed to get key %s:%s from registry!\n", section, parm);
	return NULL;
} /* get_key */


/* parses the key as string and returns a pointer */
char *get_string(reg_key *k)
{
	k->type = 0;
	return k->value;
} /* get_string */



/* parses the key as an integer and returns the value */
int get_int(reg_key *k)
{
	return (int)get_llint(k);
} /* get_int */



/* parses the key as a long integer and returns the value */
long get_long(reg_key *k)
{
	return (long)get_llint(k);
} /* get_long */


/* parses the key as a long long integer and returns the value */
long long int get_llint(reg_key *k)
{
	if (k->type != 1) {
		//k->kval.lli = strtoll(k->value, (char **) NULL, 10);      /* parse as integer */
		k->kval.lli = atoll(k->value);
		k->type = 1;
	}
	return k->kval.lli;
} /* get_llint */



/* parses the key as a float and returns the value */
float get_float(reg_key *k)
{
	return (float)get_double(k);
} /* get_float */



/* parses the key as a double and returns the value */
double get_double(reg_key *k)
{
	if (k->type != 2) {
		k->kval.dbl = atof(k->value);     /* parse as floating point */
		k->type = 2;
	}
	return k->kval.dbl;
} /* get_double */


/* writes the registry contents to a given file */
void dump_registry(registry *reg, char *filename) {
	FILE *f;
	reg_key *key;

	//fopen_s(&f, filename, "w"); // this for windows
	f = fopen(filename, "w"); // linux wants this
	if (f) {
		fprintf(f, "# Registry dump to file %s\n\n", filename);
		while (reg) {
			fprintf(f, "[%s]\n", reg->name);
			key = reg->klist;
			while (key) {
				fprintf(f, "%s = %s\n", key->name, key->value);
				key = key->next;
			} /* while */
			fprintf(f, "\n\n");
			reg = reg->next;
		} /* while */
	fclose(f);	
	} /* if */

} /* dump_registry */


int read_int_from_reg(registry *reg, char *group, char *parm)
{
	int temp;
	reg_key *k;

	k = get_key( reg, group, parm);

	if (k)
	{
		temp = get_int(k);
		return temp;
	}
	else {
		printf("UNABLE TO GET %s : %s from registry! Press any key to exit! \n", group, parm);
		getchar();
		exit(1);
	}
}


/* return true if char c is a valid token char */
int is_alpha(char c) {
	return ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || // put space to get space on paths
			(strchr("-_+*'?|!\\/()@\"$&{}<>,.;:%%", c))|| (c <= '�' && c>='�') );
} /* is_alpha */

/* end of registry.c */
