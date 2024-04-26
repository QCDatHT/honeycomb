//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//
#ifdef __cplusplus
extern "C" {
#endif

#ifndef READ_CONFIG_H
#define READ_CONFIG_H

#include <TW3EV_include/default.h>

typedef struct {
   char *content;
   uint64_t size;
   uint64_t count;
} tw3ev_string_stream_t;

typedef struct {
   tw3ev_string_stream_t **content;
   uint64_t size;
   uint64_t count;
} tw3ev_string_container_t;

tw3ev_string_stream_t *tw3ev_write_to_sstream(tw3ev_string_stream_t *restrict str, const char *restrict stuff, size_t n);
int tw3ev_write_sstream_to_file(tw3ev_string_stream_t *restrict str, const char *restrict file_path);
size_t tw3ev_write_sstream_to_file_handler(tw3ev_string_stream_t *restrict str, FILE *restrict fp);
bool tw3ev_compare_sstream_to_char(tw3ev_string_stream_t *restrict str, const char *restrict buff);
bool tw3ev_compare_sstream(tw3ev_string_stream_t *restrict str1, tw3ev_string_stream_t *restrict str2);
tw3ev_string_stream_t *tw3ev_read_file_to_sstream(const char *restrict file_path);
tw3ev_string_container_t *tw3ev_tokenize_sstream(tw3ev_string_stream_t *restrict str, const char delim[]);

void tw3ev_free_sstream(tw3ev_string_stream_t **str);
void tw3ev_clear_sstream(tw3ev_string_stream_t *str);
void tw3ev_free_scontainer(tw3ev_string_container_t **tokens);
void tw3ev_free_content_scontainer(tw3ev_string_container_t *tokens);

typedef enum { INTEGER = 0, FLOAT, EXPONENTIAL, STRING } tw3ev_input_format_e;

void tw3ev_reset_read_config();
void tw3ev_conf_append_option(const char options_local[], tw3ev_input_format_e fmt);
int tw3ev_read_configuration_file_sstream(void *(dc[]), const char fn[]);

#endif
#ifdef __cplusplus
}
#endif

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

