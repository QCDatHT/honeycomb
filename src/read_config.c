//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <TW3EV_include/read_config.h>

#define CHECK_NOT_NUM(tk, i) ((tk)[i] < '0' || (tk)[i] > '9')

tw3ev_string_stream_t *tw3ev_write_to_sstream(tw3ev_string_stream_t *  str, const char *  stuff, size_t n)
{
   if (NULL == str) {
      DA_INIT(str, tw3ev_string_stream_t, char);
   }

   size_t len = strlen(stuff);
   len = (n > len ? len : n);
   DA_APPEND_MANY(str, stuff, len);

   return str;
}

size_t tw3ev_write_sstream_to_file_handler(tw3ev_string_stream_t *  str, FILE *  fp)
{

   if (NULL == str) tw3ev_log(TW3EV_ERROR, "Trying to print NULL string stream to file.");
   if (NULL == str->content || 0 == str->size) tw3ev_log(TW3EV_ERROR, "Trying to print NULL string stream / zero-sized string stream content to file.");

   return fwrite(str->content, 1, str->count, fp);
}

int tw3ev_write_sstream_to_file(tw3ev_string_stream_t *  str, const char *  file_path)
{
   FILE *fp = fopen(file_path, "rb");
   if (!fp) tw3ev_log(TW3EV_ERROR, "Impossible to open the file <%s>", file_path);

   size_t n_char = tw3ev_write_sstream_to_file_handler(str, fp);

   fclose(fp);

   if (n_char != str->count) tw3ev_log(TW3EV_ERROR, "Wrote %ld character to file %s where count of string stream is %ld", n_char, file_path, str->count);

   return 0;
}

tw3ev_string_stream_t *tw3ev_read_file_to_sstream(const char *  file_path)
{
   FILE *fp = fopen(file_path, "rb");
   if (!fp) tw3ev_log(TW3EV_ERROR, "Impossible to open the file <%s>", file_path);

   fseek(fp, 0L, SEEK_END);
   long file_size = ftell(fp);
   rewind(fp); 

   int64_t chunks = (int64_t)file_size / ((int64_t)DA_INCREMENT);

   tw3ev_string_stream_t *sstream = (tw3ev_string_stream_t *)calloc(1, sizeof(tw3ev_string_stream_t));
   sstream->size = chunks * DA_INCREMENT + DA_INCREMENT;
   sstream->count = file_size;
   sstream->content = (char *)calloc(sstream->size, sizeof(char));

   size_t n_char = fread(sstream->content, 1, file_size, fp);
   if (n_char != (size_t)file_size) tw3ev_log(TW3EV_ERROR, "Read %ld character where size of file <%s> is %ld", n_char, file_path, file_size);
   fclose(fp);
   return sstream;
}

bool tw3ev_compare_sstream_to_char(tw3ev_string_stream_t *  str, const char *  buff)
{
   bool are_equal = true;
   size_t i = 0;
   while (are_equal) {
      are_equal = are_equal && (str->content[i] == buff[i]);
      i++;
      if (i == str->count && buff[i] != 0) are_equal = false;
      if (buff[i] == 0 && i != str->count) are_equal = false;
      if (i == str->count && buff[i] == 0) break;
   }
   return are_equal;
}

bool tw3ev_compare_sstream(tw3ev_string_stream_t *  str1, tw3ev_string_stream_t *  str2)
{
   bool are_equal = true;
   if (str1->count != str2->count) return false;

   for (size_t i = 0; i < str1->count; i++) {
      are_equal = are_equal && (str1->content[i] == str2->content[i]);
      if (!are_equal) break;
   }
   return are_equal;
}

tw3ev_string_container_t *tw3ev_tokenize_sstream(tw3ev_string_stream_t *  str, const char delim[])
{
   if (!str) tw3ev_log(TW3EV_ERROR, "Trying to tokenize null string stream");
   if (!str->content || 0 == str->size) tw3ev_log(TW3EV_ERROR, "Trying to tokenize null string stream content");
   size_t n_del = strlen(delim);
   tw3ev_string_container_t *tokens = NULL;
   DA_INIT(tokens, tw3ev_string_container_t, tw3ev_string_stream_t *);
   { 
      tw3ev_string_stream_t *temp = NULL;
      DA_INIT(temp, tw3ev_string_stream_t, char);
      DA_APPEND(tokens, temp);
   }

   bool skipping = true;

   for (size_t i = 0; i < str->count; i++) {
      bool split = false;
      for (size_t j = 0; j < n_del; j++) {
         split = split || (str->content[i] == delim[j]);
         if (split) break;
      }

      if (split && skipping) continue;
      if (split && !skipping) {
         tw3ev_string_stream_t *temp = NULL;
         DA_INIT(temp, tw3ev_string_stream_t, char);
         DA_APPEND(tokens, temp);
         skipping = true;
      }
      if (!split) {
         DA_APPEND(tokens->content[tokens->count - 1], str->content[i]);
         skipping = false;
      }
   }

   if (tokens->content[tokens->count - 1]->count == 0) {
      tw3ev_free_sstream(&(tokens->content[tokens->count - 1]));
      tokens->count--;
   }
   return tokens;
}

void tw3ev_free_sstream(tw3ev_string_stream_t **str)
{
   if (NULL == (*str)) tw3ev_log(TW3EV_ERROR, "Trying to free NULL string stream");
   if (NULL == (*str)->content && (*str)->size != 0) tw3ev_log(TW3EV_ERROR, "Trying to free NULL string stream content when size is not zero.");

   if ((*str)->size > 0) free((*str)->content);
   free((*str));
   (*str) = NULL;

   return;
}

void tw3ev_clear_sstream(tw3ev_string_stream_t *str)
{
   if (NULL == str) tw3ev_log(TW3EV_ERROR, "Trying to clear NULL string stream");
   if (NULL == str->content && str->size != 0) tw3ev_log(TW3EV_ERROR, "Trying to clear NULL string stream content when size is not zero.");

   if (str->size > 0) free(str->content);
   str->content = NULL;
   DA_INIT(str, tw3ev_string_stream_t, char);
   free(str);

   return;
}

void tw3ev_free_scontainer(tw3ev_string_container_t **tokens)
{
   if (NULL == (*tokens)) tw3ev_log(TW3EV_ERROR, "Trying to free NULL string container.");
   if (NULL == (*tokens)->content && (*tokens)->size != 0) tw3ev_log(TW3EV_ERROR, "Trying to clear NULL string container content when size is not zero.");
   if ((*tokens)->size > 0) {
      for (size_t i = 0; i < (*tokens)->count; i++)
         tw3ev_free_sstream(&((*tokens)->content[i]));
      free((*tokens)->content);
      free((*tokens));
      (*tokens) = NULL;
   }
}

void tw3ev_free_content_scontainer(tw3ev_string_container_t *tokens)
{
   if (NULL == tokens) tw3ev_log(TW3EV_ERROR, "Trying to free NULL string container.");
   if (NULL == tokens->content && 0 != tokens->size) tw3ev_log(TW3EV_ERROR, "Trying to clear NULL string container content when size is not zero.");
   if (tokens->size > 0) {
      for (size_t i = 0; i < tokens->count; i++)
         tw3ev_free_sstream(&(tokens->content[i]));
      free(tokens->content);
      tokens->content = NULL;
   }
}

#define CLEAR_LOCAL_VAR                                                                                                                                                            \
   {                                                                                                                                                                               \
      if (NULL != (options)) {                                                                                                                                                     \
         tw3ev_free_scontainer(&(options));                                                                                                                                        \
      }                                                                                                                                                                            \
      if (NULL != options_types) free(options_types);                                                                                                                              \
      options_types = NULL;                                                                                                                                                        \
      if (NULL != found_option) free(found_option);                                                                                                                                \
      found_option = NULL;                                                                                                                                                         \
   }


static tw3ev_string_container_t *options = NULL;
static tw3ev_input_format_e *options_types = NULL;
static bool *found_option = NULL;


void tw3ev_reset_read_config() { CLEAR_LOCAL_VAR; }

void tw3ev_conf_append_option(const char options_local[], tw3ev_input_format_e fmt)
{
   if (NULL == options) {
      DA_INIT(options, tw3ev_string_container_t, tw3ev_string_stream_t *);
   }

   {
      tw3ev_string_stream_t *temp = NULL;
      DA_INIT(temp, tw3ev_string_stream_t, char);
      DA_APPEND(options, temp);

      found_option = (bool *)realloc(found_option, sizeof(bool) * ((options)->count));
      options_types = (tw3ev_input_format_e *)realloc(options_types, ((options)->count) * sizeof(tw3ev_input_format_e));
   }

   size_t current = options->count - 1;

   options_types[current] = fmt;
   found_option[current] = false;
   options->content[current] = tw3ev_write_to_sstream(options->content[current], options_local, strlen(options_local));
   STR_APPEND_NULL(options->content[current]);
}

static int8_t tw3ev_check_number_nonsense(size_t opt, tw3ev_string_stream_t *tk, tw3ev_input_format_e fmt)
{
   (void)opt;
   switch (fmt) {
   case INTEGER: {
      for (size_t i = 0; i < tk->count; i++) {
         if (tk->content[i] == '\0') continue;
         if (CHECK_NOT_NUM(tk->content, i) && tk->content[0] != '+') {
            return 1;
         }
      }
      break;
   }
   case FLOAT: {

      bool dot = 0;
      if (CHECK_NOT_NUM(tk->content, 0) && tk->content[0] != '+' && tk->content[0] != '-') {
         if (tk->content[0] != '.') {
            return 1;
         } else {
            dot = true;
         }
      }
      for (size_t i = 1; i < tk->count; i++) {
         if (tk->content[i] == '\0') continue;
         if (CHECK_NOT_NUM(tk->content, i)) {
            if (tk->content[i] != '.') {
               return 1;
            } else if (!dot) {
               dot = true;
            } else {
               return 1;
            }
         }
      }
      break;
   }
   case EXPONENTIAL: {
      size_t l = tk->count;
      bool dot = false;
      if (CHECK_NOT_NUM(tk->content, 0) && tk->content[0] != '+' && tk->content[0] != '-') {
         if (tk->content[0] != '.') {
            return 1;
         } else {
            dot = true;
         }
      }
      bool e = false;
      bool sign = false;
      for (size_t i = 1; i < l; i++) {
         if (tk->content[i] == '\0') continue;
         if (CHECK_NOT_NUM(tk->content, i)) {
            if (tk->content[i] != '.' && tk->content[i] != 'e' && tk->content[i] != 'E' && tk->content[i] != '+' && tk->content[i] != '-') {
               return 1;
            } else {
               switch (tk->content[i]) {
               case '.': {
                  if (!dot) {
                     dot = true;
                     break;
                  } else goto error_goto;
                  break;
               }
               case '+': {
                  if (!sign) {
                     sign = true;
                     break;
                  } else goto error_goto;
                  break;
               }
               case '-': {
                  if (!sign) {
                     sign = true;
                     break;
                  } else goto error_goto;
                  break;
               }
               case 'e':
               case 'E': {
                  if (!e) {
                     e = true;
                     break;
                  } else goto error_goto;
                  break;
               }

               default: {
               error_goto:
                  return 1;
               }
               }
            }
         }
      }
      break;
   }
   default:
      break;
   }
   return 0;
}

static void tw3ev_fill_variable(size_t opt, tw3ev_string_stream_t *tk, void *value)
{
   if (opt >= options->count) tw3ev_log(TW3EV_ERROR, "Out of bound option.");

   int8_t is_not_ok = 0;
   if (options_types[opt] == FLOAT || options_types[opt] == EXPONENTIAL) {
      is_not_ok = tw3ev_check_number_nonsense(opt, tk, FLOAT);
      is_not_ok = is_not_ok && tw3ev_check_number_nonsense(opt, tk, EXPONENTIAL);
      if (is_not_ok) {
         tw3ev_log(TW3EV_INFO,
                   "Something wrong in the format for %s. Expected a "
                   "floating point number, instead found <<%s>>",
                   options->content[opt]->content, tk->content);
         tw3ev_log(TW3EV_INFO, "Resorting to default value, see `config.out\'");
         return;
      }

   } else if (options_types[opt] == INTEGER) {
      is_not_ok = tw3ev_check_number_nonsense(opt, tk, INTEGER);
      if (is_not_ok) {
         tw3ev_log(TW3EV_INFO,
                   "Something wrong in the format for %s. Expected an "
                   "integer number, instead found <<%s>>",
                   options->content[opt]->content, tk->content);
         tw3ev_log(TW3EV_INFO, "Resorting to default value, see `config.out\'");
         return;
      }
   }

   found_option[opt] = true;

   switch (options_types[opt]) {
   case INTEGER: {
      int32_t *location = (int32_t *)value;
      if (NULL == location) tw3ev_log(TW3EV_ERROR, "Trying to write integer to NULL location.");
      *location = atoi(tk->content);
      break;
   }
   case FLOAT: {
      double *location = (double *)value;
      if (NULL == location) tw3ev_log(TW3EV_ERROR, "Trying to write float to NULL location.");
      *location = strtod(tk->content, NULL);
      break;
   }
   case EXPONENTIAL: {
      double *location = (double *)value;
      if (NULL == location) tw3ev_log(TW3EV_ERROR, "Trying to write exponential float to NULL location.");
      *location = strtod(tk->content, NULL);
      break;
   }
   case STRING: {
      char *temp = (char *)value;
      if (NULL == temp) tw3ev_log(TW3EV_ERROR, "Trying to write string to NULL location.");
      size_t how_many = min2((uint64_t)(BASEFOLDER_MAX_L + 1), tk->count);
      if (how_many != tk->count) tw3ev_log(TW3EV_WARNING, "Cutting token to fit in standard string size. Please check if this is intended behavior.");
      for (size_t i = 0; i < how_many; i++)
         temp[i] = tk->content[i];

      break;
   }
   default:
      break;
   }

   return;
}

static int8_t tw3ev_check_for_skipping_in_sstream(tw3ev_string_stream_t *str)
{
   for (size_t i = 0; i < str->count; i++) {
      if (str->content[i] == '\n' || str->content[i] == '#' || str->content[i] == '\0') return 2;
      if (str->content[i] == ' ' || str->content[i] == '\t') {
         (str->content)++;
         str->count--;
         i++;
      } else {
         return 0;
      }
   }
   return 1;
}

int tw3ev_read_configuration_file_sstream(void *(dc[]), const char fn[])
{
   tw3ev_string_stream_t *full_file = tw3ev_read_file_to_sstream(fn);
   tw3ev_string_container_t *lines = tw3ev_tokenize_sstream(full_file, "\n");
   tw3ev_free_sstream(&full_file);
   char delim[24] = " ,;:\n\t=";

   for (size_t i = 0; i < lines->count; i++) {
      tw3ev_string_stream_t *line = lines->content[i];
      if (tw3ev_check_for_skipping_in_sstream(line) != 0) continue;

      for (size_t j = 0; j < line->count; j++) 
         if (line->content[j] == '#') line->count = j;
      

      tw3ev_string_container_t *tokens = tw3ev_tokenize_sstream(line, delim);
      STR_APPEND_NULL(tokens->content[0]);
      if (tokens->count <= 1) {
         tw3ev_log(TW3EV_WARNING, "In reading the configuration file `%s\': empty token value <<%s>>", fn, tokens->content[0]->content);
         tw3ev_free_scontainer(&tokens);
         continue;
      }

      size_t opt = 0;
      for (opt = 0; opt < options->count; opt++) {
         if (tw3ev_compare_sstream(tokens->content[0], options->content[opt])) break;
      }

      if (opt >= options->count) {
         tw3ev_log(TW3EV_WARNING, "In reading the configuration file `%s\': unrecognize token <<%s>>", fn, tokens->content[0]->content);
         tw3ev_free_scontainer(&tokens);
         continue;
      }

      STR_APPEND_NULL(tokens->content[1]);

      tw3ev_fill_variable(opt, tokens->content[1], dc[opt]);

      for (size_t j = 2; j < tokens->count; j++) {
         STR_APPEND_NULL(tokens->content[j]);
         tw3ev_log(TW3EV_WARNING, "Ignored tokens %s", tokens->content[j]->content);
      }
      tw3ev_free_scontainer(&tokens);
   }
   for (size_t i = 0; i < options->count; i++) {
      if (!found_option[i]) {
         STR_APPEND_NULL(options->content[i]);
         tw3ev_log(TW3EV_WARNING, "Option not found or incorrect format: %s", options->content[i]->content);
      }
   }

   return 0;
}


 // Copyright (C) 2024 Simone Rodini; Lorenzo Rossi
 // This program is free software; you can redistribute it and/or modify
 // it under the terms of the GNU General Public License as published by
 // the Free Software Foundation; either version 2 of the License, or
 // (at your option) any later version.
 // 
 // This program is distributed in the hope that it will be useful,
 // but WITHOUT ANY WARRANTY; without even the implied warranty of
 // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 // GNU General Public License for more details.
 // 
 // You should have received a copy of the GNU General Public License along
 // with this program; if not, write to the Free Software Foundation, Inc.,
 // 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

