#include <stdio.h>
#include <stdlib.h>
#include "ezxml/ezxml.h"


/* 
 * Interface routines for building streams at run-time; defined in mpas_stream_manager.F
 */
void stream_mgr_init_c(char *, int *);


void xml_stream_parser(char *fname, void *manager)
{
	FILE *iofd;
	ezxml_t streams;
	ezxml_t foo_xml;
	const char *attr;
	int err;

	fprintf(stderr, "MGD DEV begin parsing run-time I/O from %s\n", fname);

	iofd = fopen(fname, "r");
	if (!iofd) {
		fprintf(stderr, "Error: Could not open run-time I/O config file %s\n", fname);
		return;
	}

	streams = ezxml_parse_fp(iofd);
	if (!streams) {
		fprintf(stderr, "Error: Problems encountered while parsing run-time I/O config file %s\n", fname);
		return;
	}	

	err = 0;
	for (foo_xml = ezxml_child(streams, "stream"); foo_xml; foo_xml = ezxml_next(foo_xml)) {
		attr = ezxml_attr(foo_xml, "name");
		fprintf(stderr, "MGD DEV Found stream named %s\n", attr);
		stream_mgr_create_stream_c(manager, attr, &err);
		err++;
	}

	fprintf(stderr, "MGD DEV done parsing run-time I/O\n");
}
