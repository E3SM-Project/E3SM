#include <stdio.h>
#include <stdlib.h>
#include <regex.h>

#define MAX_LEN 1024

void check_regex_match(const char * pattern, const char * str, int *imatch){
	regex_t regex;
	char bracketed_pattern[MAX_LEN];
	int ierr, len;

	*imatch = 0;
	len = snprintf(bracketed_pattern, 1024, "^%s$", pattern);
	if ( len >= MAX_LEN ) {
		*imatch = -1;
		return;
	}

	ierr = regcomp(&regex, bracketed_pattern, 0);
	if ( ierr ) {
		*imatch = -1;
		return;
	}

	ierr = regexec(&regex, str, 0, NULL, 0);

	regfree(&regex);

	if ( !ierr ) {
		*imatch = 1;
	} else if ( ierr == REG_NOMATCH ) {
		*imatch = 0;
	} else {
		*imatch = -1;
	}
}

