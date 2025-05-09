//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <honeycomb/integration.h>

// NOTE: default to GK21 if nothing is specified
static integration_rule *local_integration_rule = integration_rule_gk61;

void set_integration_routine(integration_rule *ir) { local_integration_rule = ir; }

void integration(double (*fnc)(double, void *), void *p_fnc, double a, double b, double epsabs, double *result) { local_integration_rule(fnc, p_fnc, a, b, epsabs, result); }

void integration_rule_gk21(double (*fnc)(double, void *), void *p_fnc, double a, double b, double epsabs, double *result)
{
   integration_qag(fnc, p_fnc, a, b, epsabs, result, RULE_TW3EV_GK21);
}

void integration_rule_gk31(double (*fnc)(double, void *), void *p_fnc, double a, double b, double epsabs, double *result)
{
   integration_qag(fnc, p_fnc, a, b, epsabs, result, RULE_TW3EV_GK31);
}

void integration_rule_gk41(double (*fnc)(double, void *), void *p_fnc, double a, double b, double epsabs, double *result)
{
   integration_qag(fnc, p_fnc, a, b, epsabs, result, RULE_TW3EV_GK41);
}

void integration_rule_gk61(double (*fnc)(double, void *), void *p_fnc, double a, double b, double epsabs, double *result)
{
   integration_qag(fnc, p_fnc, a, b, epsabs, result, RULE_TW3EV_GK61);
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

