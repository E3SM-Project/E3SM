// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//
#include <stdio.h>
#include <stdarg.h>

#ifdef FORTPRINTF_TESTING
#define MAX_LINE_LEN 20
#else
#define MAX_LINE_LEN 132
#endif

char printbuf[MAX_LINE_LEN+2];
char fbuffer[1024];
int nbuf = 0;

/*
 * TODO: Split quoted strings wherever necessary to fit within line length limit, rather
 *       than worrying about splitting strings at spaces
 */

int fortprintf(FILE * fd, char * str, ...)/*{{{*/
{
	int i, nl, sp, sp_inquotes, inquotes, q;
	int lastchar;
	int errorcode;
	va_list ap;

#ifdef FORTPRINTF_DEBUG
	printf("call to fortprintf\n");
	printf("Format string is: %s\n", str);
	printf("\n");
#endif

	/* Assume no errors */
	errorcode = 0;

	/* Add formatted string to the buffer of fortran code to be written */
	va_start(ap, str);
	i = vsnprintf(fbuffer+nbuf, 1024-nbuf, str, ap);
#ifdef FORTPRINTF_DEBUG
	printf("Full buffer is:\n");
	printf("%s\n", fbuffer);
#endif
	va_end(ap);

	/* Set the next free position in the fortran buffer */
	nbuf = nbuf + i;

	inquotes = 0;
	q = -1;

	do {

		nl = sp = -1;

		/* Scan through the max line length - 1 (since we may have to add an & character) or the end of the buffer, whichever comes first */
		for (i=0; i<MAX_LINE_LEN-1 && i<nbuf; i++) {
			lastchar = (i == nbuf-1) ? 1 : 0;
			if (fbuffer[i] == '\'') {
				if ( fbuffer[i+1] != '\'' || lastchar ) {
					inquotes = (inquotes + 1) % 2;
					q = inquotes ? i : -1;
				} else if ( !lastchar && fbuffer[i+1] == '\'' ) {
					i++;
				}
			}
			if (fbuffer[i] == '\n') nl = i;                                                               /* The last occurrence of a newline */
			if (fbuffer[i] == ' ' && !lastchar && fbuffer[i+1] != '&') {                             /* The last occurrence of a space */
				sp = i;
				sp_inquotes = inquotes;
			}
		}

#ifdef FORTPRINTF_DEBUG
		printf("1: i = %d, nl = %d, sp = %d, fbuffer[i] = %c, inquotes = %d\n", i, nl, sp, fbuffer[i], inquotes);
#endif

		/* If we haven't reached the column limit, don't consider breaking the line yet */
		if (nbuf <= MAX_LINE_LEN) sp = -1;

#ifdef FORTPRINTF_DEBUG
		printf("2: i = %d, nl = %d, sp = %d, fbuffer[i] = %c\n", i, nl, sp, fbuffer[i]);
#endif

		/* If the charater at column MAX_LINE_LEN happens to be a newline or the only space, though, we mark it */
		if (i == MAX_LINE_LEN-1) {
			if (fbuffer[i] == '\n') nl = i;
			if (fbuffer[i] == ' ' && sp < 0)  sp = i-1;   /* Insert & before the space so it will fit within column limit */
		}

#ifdef FORTPRINTF_DEBUG
		printf("3: i = %d, nl = %d, sp = %d, fbuffer[i] = %c\n", i, nl, sp, fbuffer[i]);
#endif

		/* If we have a newline */
		if (nl >= 0) {
			snprintf(printbuf, nl+2, "%s", fbuffer);
			fprintf(fd, "%s", printbuf);
			nl++;

			/* Shift unprinted contents of fortran buffer to the beginning */
			for (i=0; nl<nbuf; i++, nl++)
				fbuffer[i] = fbuffer[nl];
			nbuf = i;
		}
		/* Else if we found a place to break the line */
		else if (sp >= 0) {
			snprintf(printbuf, sp+2, "%s", fbuffer);
			i = sp+1;
			if (sp_inquotes && (sp > q)) printbuf[i++] = '\'';
			printbuf[i++] = '&';
			printbuf[i++] = '\n';
			printbuf[i++] = '\0';
			fprintf(fd, "%s", printbuf);
			sp++;
			i = 0;
			if (sp_inquotes && (sp > q)) {
				inquotes = (inquotes + 1) % 2;
				fbuffer[i++] = '/';
				fbuffer[i++] = '/';
				fbuffer[i++] = '\'';
			}

			/* Shift unprinted contents of fortran buffer to the beginning */
			for ( ; sp<nbuf; i++, sp++)
				fbuffer[i] = fbuffer[sp];
			nbuf = i;
		}
		else if (nbuf > MAX_LINE_LEN) {
			fprintf(fd, "!!! Error: Output could not be formatted by fortprintf() !!!\n");
			nbuf = 0;
			errorcode = 1;
		}

	} while (nl >= 0 || sp >= 0);

	return errorcode;
}/*}}}*/

void fortprint_flush(FILE * fd)/*{{{*/
{
	snprintf(printbuf, nbuf+1, "%s", fbuffer);
	fprintf(fd, "%s", printbuf);
	nbuf = 0;
}/*}}}*/

#ifdef FORTPRINTF_TESTING
void print_result(int test_num, int err_code)/*{{{*/
{
	if(err_code == 0){
		printf("Test %02d:   PASSED\n", test_num);
	} else {
		printf("Test %02d: **FAILED\n", test_num);
	}
}/*}}}*/

int main()/*{{{*/
{
	FILE * foo;
	int err;

	/* Tests writing a line with NO space that is below the column limit */
	foo = fopen("test01.inc","w");
	err = fortprintf(foo, "123456789\n");
	print_result(1, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests writing a line with space that is below the column limit */
	foo = fopen("test02.inc","w");
	err = fortprintf(foo, "12345 789\n");
	print_result(2, err);
	fortprint_flush(foo);
	fclose(foo);

	/*** Test lines that are less than 20 chars long ***/

	/* Tests the case where we write a newline at the last column with chances to break the line earlier */
	foo = fopen("test03.inc","w");
	err = fortprintf(foo, "123456789 12345678\n");
	print_result(3, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests the case where we write a newline at the last column with NO chances to break the line earlier */
	foo = fopen("test04.inc","w");
	err = fortprintf(foo, "123456789012345678\n");
	print_result(4, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with a space occurring in the NEXT-TO-LAST column plus another line */
	foo = fopen("test05.inc","w");
	err = fortprintf(foo, "123456789 1234567 0");
	err += fortprintf(foo, "1234\n");
	print_result(5, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with a space occurring in the LAST column plus another line */
	foo = fopen("test06.inc","w");
	err = fortprintf(foo, "123456789 12345678 ");
	err += fortprintf(foo, "1234\n");
	print_result(6, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with the first space occurring in the NEXT-TO-LAST column plus another line */
	foo = fopen("test07.inc","w");
	err = fortprintf(foo, "12345678901234567 0");
	err += fortprintf(foo, "1234\n");
	print_result(7, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with the first space occurring in the LAST column plus another line */
	foo = fopen("test08.inc","w");
	err = fortprintf(foo, "123456789012345678 ");
	err += fortprintf(foo, "1234\n");
	print_result(8, err);
	fortprint_flush(foo);
	fclose(foo);

	/*** Test lines that are exactly 20 chars long ***/

	/* Tests the case where we write a newline at the last column with chances to break the line earlier */
	foo = fopen("test09.inc","w");
	err = fortprintf(foo, "123456789 123456789\n");
	print_result(9, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests the case where we write a newline at the last column with NO chances to break the line earlier */
	foo = fopen("test10.inc","w");
	err = fortprintf(foo, "1234567890123456789\n");
	print_result(10, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with a space occurring in the NEXT-TO-LAST column plus another line */
	foo = fopen("test11.inc","w");
	err = fortprintf(foo, "123456789 12345678 0");
	err += fortprintf(foo, "1234\n");
	print_result(11, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with a space occurring in the LAST column plus another line */
	foo = fopen("test12.inc","w");
	err = fortprintf(foo, "123456789 123456780 ");
	err += fortprintf(foo, "1234\n");
	print_result(12, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with the first space occurring in the NEXT-TO-LAST column plus another line */
	foo = fopen("test13.inc","w");
	err = fortprintf(foo, "123456789012345678 0");
	err += fortprintf(foo, "1234\n");
	print_result(13, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with the first space occurring in the LAST column plus another line */
	foo = fopen("test14.inc","w");
	err = fortprintf(foo, "1234567890123456780 ");
	err += fortprintf(foo, "1234\n");
	print_result(14, err);
	fortprint_flush(foo);
	fclose(foo);

	/*** Test lines that are more than 21 chars long ***/

	/* Tests the case where we write a newline at the last column with chances to break the line earlier */
	foo = fopen("test15.inc","w");
	err = fortprintf(foo, "123456789 1234567890\n");
	print_result(15, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests the case where we write a newline at the last column with NO chances to break the line earlier */
	foo = fopen("test16.inc","w");
	err = fortprintf(foo, "1234567890123456789\n");
	print_result(16, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with a space occurring in the NEXT-TO-LAST column plus another line */
	foo = fopen("test17.inc","w");
	err = fortprintf(foo, "123456789 123456789 0");
	err = fortprintf(foo, "1234\n");
	print_result(17, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with a space occurring in the LAST column plus another line */
	foo = fopen("test18.inc","w");
	err = fortprintf(foo, "123456789 1234567890 ");
	err += fortprintf(foo, "1234\n");
	print_result(18, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with the first space occurring in the NEXT-TO-LAST column plus another line */
	foo = fopen("test19.inc","w");
	err = fortprintf(foo, "1234567890123456789 0");
	err += fortprintf(foo, "1234\n");
	print_result(19, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with the first space occurring in the LAST column plus another line */
	foo = fopen("test20.inc","w");
	err = fortprintf(foo, "12345678901234567890 ");
	err += fortprintf(foo, "1234\n");
	print_result(20, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with single quotes */
	foo = fopen("test21.inc", "w");
	err = fortprintf(foo, "\'1234567890\'");
	print_result(21, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with double quotes */
	foo = fopen("test22.inc", "w");
	err = fortprintf(foo, "\'1234\'\'56789\'");
	print_result(22, err);
	fortprint_flush(foo);
	fclose(foo);

	/* Tests a line with double quotes and a line break after the double quotes */
	foo = fopen("test22.inc", "w");
	err = fortprintf(foo, "\'1234567\'\'890 1234567890\'");
	print_result(22, err);
	fortprint_flush(foo);
	fclose(foo);


	return 0;
}/*}}}*/
#endif
