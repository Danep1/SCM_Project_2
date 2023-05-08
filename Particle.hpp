#pragma once
#include <random>

#include "vector.hpp"

number_t distance(const r_point & v1, const r_point & v2);

class Particle
{
private:
	number_t m_mass;
	number_t m_charge;
	r_point m_pos;
	r_point m_prev_pos;
	r_point m_v;

	number_t m_viscosity = 1.0;
	number_t m_random_function_factor = 0.0; // 

public:
	r_point m_a;

public:
	Particle() noexcept = default;

	explicit Particle(number_t mass, number_t charge, const r_point& pos, const r_point& prev_pos, const r_point& v) noexcept :
		m_mass(mass), m_charge(charge), m_pos(pos), m_v(v), m_prev_pos(prev_pos), m_a(r_point(0.0, 0.0, 0.0))
	{}

	~Particle() noexcept = default;

public:
	const auto get_pos() inline const noexcept
	{
		return m_pos;
	}

	const auto get_prev_pos() inline const noexcept
	{
		return m_prev_pos;
	}

	const auto get_v() inline const noexcept
	{
		return m_v;
	}

	const auto get_mass() inline const noexcept
	{
		return m_mass;
	}

	const auto get_charge() inline const noexcept
	{
		return m_charge;
	}

	const auto get_T() inline const noexcept
	{
		return m_mass * (m_v * m_v) / 2.0 * 1.036427061146e-11; // eV
	}

	void save_cur_pos_as_prev() 
	{
		m_prev_pos = m_pos;
	}

	void update(number_t dt) noexcept;

	void move_with(const r_point& dr) noexcept;
};