#pragma once
#include "Cell.hpp"
#include <fstream>
#include <iostream>

inline const number_t Kb = 8.617333262e-5;

class Processor
{
private:

	const number_t m_const_a = 4.1; // ангстрем

	const number_t m_length = m_const_a * 5U;	// ангстрем
	const number_t m_width = m_const_a * 5U;	// ангстрем
	const number_t m_height = m_const_a * 15U;	// ангстрем
	
	const number_t m_v_max = 1.0e+0;

	const std::size_t m_N_particls_in_row = 5U;

	const number_t m_sigma = 2.6175; // m_const_a / 1.095 / std::sqrt(2.0); // ангстрем
	const number_t m_R_cut = m_sigma * 2.5;   // ангстрем
	const number_t m_U_0 = 4577.6 * Kb; // eV

	const number_t t = 1.0e+1;	// pikosec
	const number_t dt = 1.0e+0; // pikosec
	const std::size_t m_N_steps = t / dt;
	const std::size_t m_N_update = 1U;

	Cell m_cell;

public:
	Processor() : m_cell(m_sigma, m_R_cut, m_U_0, m_v_max, r_point(m_length, m_width, m_height), m_N_particls_in_row, dt) {}

	~Processor() noexcept = default;

	void start();

	void write_energy(std::ofstream& fstream, number_t t);

	void write_current_system(std::ofstream& fstream, number_t t);
};