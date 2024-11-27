/*  read_fast.c --  showcases the htslib api usage

    Copyright (C) 2023 Genome Research Ltd.

    Author: Vasudeva Sarma <vasudeva.sarma@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE

*/

/* The purpose of this code is to demonstrate the library apis and need proper error handling and optimisation */

#include <getopt.h>
#include <unistd.h>
#include <htslib/sam.h>

/// print_usage - show usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void print_usage(FILE *fp)
{
    fprintf(fp, "Usage: read_fast <infile>\n\
Reads the fasta/fastq file and shows the content.\n");
    return;
}

/// main_demo - start of the demo
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main(int argc, char *argv[])
{
    const char *inname = NULL;              //input file name
    int c = 0, ret = EXIT_FAILURE;
    samFile *infile = NULL;                 //sam file
    sam_hdr_t *in_samhdr = NULL;            //header of file
    bam1_t *bamdata = NULL;                 //to hold the read data

    if (argc != 2) {
        print_usage(stdout);
        goto end;
    }
    inname = argv[1];

    //initialize
    if (!(bamdata = bam_init1())) {
        printf("Failed to initialize bamdata\n");
        goto end;
    }
    //open input files - r reading
    if (!(infile = sam_open(inname, "r"))) {
        printf("Could not open %s\n", inname);
        goto end;
    }
    if (infile->format.format != fasta_format && infile->format.format != fastq_format) {
        printf("Invalid file specified\n");
        goto end;
    }

    //read header
    if (!(in_samhdr = sam_hdr_read(infile))) {
        printf( "Failed to read header from file\n");
        goto end;
    }

    //read data
    while ((c = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        printf("\nname: ");
        printf("%s", bam_get_qname(bamdata));
        printf("\nsequence: ");
        for (c = 0; c < bamdata->core.l_qseq; ++c) {
            printf("%c", seq_nt16_str[bam_seqi(bam_get_seq(bamdata), c)]);
        }
        if (infile->format.format == fastq_format) {
            printf("\nquality: ");
            for (c = 0; c < bamdata->core.l_qseq; ++c) {
                printf("%c", bam_get_qual(bamdata)[c] + 33);
            }
        }
    }
    printf("\n");
    if (c != -1) {
        //error
        printf("Failed to get data\n");
        goto end;
    }
    //else -1 / EOF
    ret = EXIT_SUCCESS;
end:
    //clean up
    if (in_samhdr) {
        sam_hdr_destroy(in_samhdr);
    }
    if (infile) {
        sam_close(infile);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }
    return ret;
}