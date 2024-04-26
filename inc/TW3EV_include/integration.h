//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//
#ifdef __cplusplus
extern "C" {
#endif

#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <TW3EV_include/default.h>

typedef enum { RULE_TW3EV_GK21, RULE_TW3EV_GK31, RULE_TW3EV_GK41, RULE_TW3EV_GK61 } tw3ev_integration_rule_e;

typedef void integration_rule(double (*fnc)(double, void *), void *p_fnc, double a, double b, double epsabs, double *result);
void integration_qag(double (*fnc)(double, void *), void *p_fnc, double a, double b, double epsabs, double *result, tw3ev_integration_rule_e gk_type);

void integration(double (*fnc)(double, void *), void *p_fnc, double a, double b, double epsabs, double *result);
void set_integration_routine(integration_rule *ir);

integration_rule integration_rule_gk21;
integration_rule integration_rule_gk31;
integration_rule integration_rule_gk41;
integration_rule integration_rule_gk61;

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

