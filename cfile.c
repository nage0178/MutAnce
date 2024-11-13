/*
    Copyright (C) 2016-2022 Tomas Flouri, Bruce Rannala and Ziheng Yang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <t.flouris@ucl.ac.uk>,
    Department of Genetics, Evolution and Environment,
    University College London, Gower Street, London WC1E 6BT, England
*/

#include "MutAnce.h"

static char buffer[LINEALLOC];
static char * line = NULL;
static size_t line_size = 0;
static size_t line_maxsize = 0;


static void reallocline(size_t newmaxsize)
{
  char * temp = (char *)xmalloc((size_t)newmaxsize*sizeof(char));

  if (line)
  {
    memcpy(temp,line,line_size*sizeof(char));
    free(line);
  }
  line = temp;
  line_maxsize = newmaxsize;
}

static char * getnextline(FILE * fp)
{
  size_t len = 0;

  line_size = 0;

  /* read from file until newline or eof */
  while (fgets(buffer, LINEALLOC, fp))
  {
    len = strlen(buffer);

    if (line_size + len > line_maxsize)
      reallocline(line_maxsize + LINEALLOC);

    memcpy(line+line_size,buffer,len*sizeof(char));
    line_size += len;

    if (buffer[len-1] == '\n')
    {
        line[line_size-1] = 0;

      return line;
    }
  }

  if (!line_size)
  {
    free(line);
    line_maxsize = 0;
    line = NULL;
    return NULL;
  }

  if (line_size == line_maxsize)
    reallocline(line_maxsize+1);

  line[line_size] = 0;
  return line;
}

static long is_emptyline(const char * line)
{
  size_t ws = strspn(line, " \t\r\n");
  if (!line[ws] || line[ws] == '*' || line[ws] == '#') return 1;
  return 0;
}


static long get_string(const char * line, char ** value)
{
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except star and hash */
  while (*p && *p != '*' && *p != '#') p++;


  /* now go back, skipping all white-space until a character occurs */
  for (--p; *p == ' ' || *p == '\t' || *p == '\r' || *p =='\n'; --p);
  
  char * end = p+1;
  *end = 0;

  *value = xstrdup(start);
  free(s);

  return ws + end - start;
}

static long get_long(const char * line, long * value)
{
  int ret,len=0;
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except star, hash and whitespace */
  char * end = start + strcspn(start," \t\r\n*#");

  *end = 0;

  ret = sscanf(start, "%ld%n", value, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(start)))
  {
    free(s);
    return 0;
  }

  free(s);
  return ws + end - start;
}



static long parse_long(const char * line, long * value)
{
  long ret = 0;
  char * s = xstrdup(line);
  char * p = s;

  long count;

  count = get_long(p, value);
  if (!count) goto l_unwind;

  p += count;

  if (is_emptyline(p)) ret = 1;

l_unwind:
  free(s);
  return ret;
}

static long get_token(char * line, char ** token, char ** value)
{
  char * p = line;

  /* here we parse lines which may be in the form:
     
     token = value
  */

  /* skip all white-space */
  while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') ++p;

  /* is it a blank line or comment ? */
  if (!*p || *p == '*' || *p == '#') return 0;

  /* store address of token's beginning */
  *token = p;

  /* find occurrence of '=' */
  while (*p && *p != '=') ++p;

  /* if no '=' found return error */
  if (!*p) return -1;

  /* '=' was found, store pointer to value */
  *value = p+1;

  /* '=' was found, move back to the last letter of token (ignore whitespace) */
  for (--p; *p == ' ' || *p == '\t' || *p == '\r' || *p =='\n'; --p);

  /* return length of token */
  return p - *token + 1;
}


static void check_validity(char * opt_outfile, char * opt_msafile, char * opt_gfile, char * opt_mfile,
		long locusSim)
{
  /*if (!opt_streenewick)
    fatal("Initial species tree newick format is required in 'species&tree'"); */

/*if (!opt_outfile)
    fatal("Option 'outfile' is required"); */

  if (!opt_msafile)
    fatal("Option 'seqfile' is required");

  if (!opt_gfile)
    fatal("Option 'treefile' is required");

  if (!opt_mfile)
    fatal("Option 'migfile' is required");

  if (locusSim < 1) {
    fatal("Option 'locus' must be a positive integer greater than zero");
  }
  
/*  if (opt_samples < 1)
    fatal("Option 'nsample' must be a positive integer greater than zero");
  
  if (opt_samplefreq < 1)
    fatal("Option 'sampfreq' must be a positive integer greater than zero");
*/

}


void load_cfile(char * opt_cfile, char ** opt_msafile, 
		char ** opt_outfile, 
		char ** opt_gfile, char ** opt_mfile,
		char ** opt_subfile,
		long * opt_seed, long * locusSim, long *opt_site)
{
  long line_count = 0;
  FILE * fp;

  /* the following variable is used for checking whether we have a newick
     string in the species&tree tag, in the case of 1 species. For species
     trees we do not accept a tree, whereas for network we require a newick
     string. The program always reads a line. If that line is a tree it is
     processed, otherwise this variable is set such that we do not read another
     line */
  long line_not_processed = 0;

  fp = xopen(opt_cfile,"r");

  while (line_not_processed || getnextline(fp))
  {
    int valid = 0;
    char * token;
    char * value;
    long token_len;

    line_not_processed = 0;

    ++line_count;
    token_len = get_token(line,&token,&value);

    if (!token_len) continue;
    if (token_len < 0)
      fatal("Invalid syntax when parsing file %s on line %ld",
            opt_cfile, line_count);
    
    if (token_len == 4)
    {
      if (!strncasecmp(token,"seed",4))
      {
        if (!parse_long(value,opt_seed))
          fatal("Option 'seed' expects one integer (line %ld)", line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"site",4)) {
        if (!parse_long(value,opt_site) || opt_site < 0)
          fatal("Option 'site' expects one positive integer (line %ld)", line_count);
        valid = 1;
      }
    }
    else if (token_len == 5)
    {
      if (!strncasecmp(token,"nloci",5))
      {
        if (!parse_long(value,&opt_locus_count) || opt_locus_count < 0)
          fatal("Option 'nloci' expects a positive integer or zero (line %ld)",
                line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"locus",5))
      {
        if (!get_long(value,locusSim))
          fatal("Option %s expects a long (line %ld)", token, line_count);
        valid = 1;
      }
    }
    else if (token_len == 7)
    {
      if (!strncasecmp(token,"seqfile",7))
      {
        if (!get_string(value, opt_msafile))
          fatal("Option %s expects a string (line %ld)", token, line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"outfile",7))
      {
        if (!get_string(value, opt_outfile))
          fatal("Option %s expects a string (line %ld)", token, line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"migfile",7))
      {
        if (!get_string(value, opt_mfile))
          fatal("Option %s expects a string (line %ld)", token, line_count);
        valid = 1;
      }
      else if (!strncasecmp(token,"HKYfile",7))
      {
        if (!get_string(value, opt_subfile))
          fatal("Option %s expects a string (line %ld)", token, line_count);
        valid = 1;
      }
    } else if (token_len == 9) {

      if (!strncasecmp(token,"gtreefile",8))
      {
        if (!get_string(value,opt_gfile))
          fatal("Option %s expects a string (line %ld)", token, line_count);
        valid = 1;
      }
    }
    if (!valid)
      fatal("Invalid syntax when parsing file %s on line %ld",
            opt_cfile, line_count);
  }

  fclose(fp);
  check_validity(*opt_outfile, *opt_msafile, *opt_gfile, *opt_mfile, *locusSim );

}

//int parsefile_doubles(const char * filename,
//                      long n,
//                      double * outbuffer,
//                      long * errcontext)
//{
//  long line_count = 0;
//  long count;
//  FILE * fp = xopen(filename,"r");
//
//  long entry = 0;
//
//  while (getnextline(fp))
//  {
//    ++line_count;
//
//    char * p = line;
//
//    while (!is_emptyline(p))
//    {
//      if (entry == n)
//      {
//        fclose(fp);
//        return ERROR_PARSE_MORETHANEXPECTED;
//      }
//      count = get_double(p, outbuffer+entry++);
//      if (!count)
//      {
//        fclose(fp);
//        *errcontext = line_count;
//        return ERROR_PARSE_INCORRECTFORMAT;
//      }
//
//      p += count;
//      
//    }
//  }
//  if (entry != n)
//  {
//    fclose(fp);
//    *errcontext = entry;
//    return ERROR_PARSE_LESSTHANEXPECTED;
//  }
//  
//  fclose(fp);
//  return 0;
//}
