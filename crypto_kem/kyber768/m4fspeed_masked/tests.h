/* Copyright 2022 UCLouvain, Belgium and PQM4 contributors
 *
 * This file is part of pqm4_masked.
 *
 * pqm4_masked is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, version 3.
 *
 * pqm4_masked is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * pqm4_masked. If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef TESTS_H
#define TESTS_H
unsigned int test_function();
unsigned int test_convertions_APoly();
unsigned int test_convertions_bitslice();
unsigned int test_xor_bitslice();
unsigned int test_and_bitslice();
unsigned int test_secadd();
unsigned int test_secadd_modp();
unsigned int test_secadd_constant_bmsk();
unsigned int test_secadd_constant();
unsigned int test_seca2b();
unsigned int test_seca2b_modp();
unsigned int test_secb2a_modp();
unsigned int test_secb2a_modp_1bit();
unsigned int test_seccompress();
unsigned int test_cbd();
#endif
