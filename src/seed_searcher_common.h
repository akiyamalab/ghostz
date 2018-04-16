/*
 * seed_searcher_common.h
 *
 *  Created on: 2014/05/13
 *      Author: shu
 */

#ifndef SEED_SEARCHER_COMMON_H_
#define SEED_SEARCHER_COMMON_H_

#include "reduced_alphabet_variable_hash_function.h"

class SeedSearcherCommon {
public:
	typedef int Distance;
	typedef ReducedAlphabetVariableHashFunction HashFunction;
	typedef HashFunction::Hash Hash;

	static SeedSearcherCommon::Distance CalculateDistance(
			const AlphabetCoder::Code *subsequence0_center,
			const AlphabetCoder::Code *subsequence1_center,
			const uint32_t sequence_length,
			const std::vector<AlphabetCoder::Code> &redueced_code_map);

	static uint32_t CalculateSeedPosition(uint32_t hashed_sequence_start,
			uint32_t hashed_sequence_length);

private:
};

inline SeedSearcherCommon::Distance SeedSearcherCommon::CalculateDistance(
		const AlphabetCoder::Code *subsequence0_center,
		const AlphabetCoder::Code *subsequence1_center,
		const uint32_t sequence_length,
		const std::vector<AlphabetCoder::Code> &redueced_code_map) {
	const uint32_t foward_direction = sequence_length / 2;
	uint32_t mismatch_count = 0;
	const AlphabetCoder::Code *subsequence0 = subsequence0_center
			- foward_direction;
	const AlphabetCoder::Code *subsequence1 = subsequence1_center
			- foward_direction;
	for (uint32_t i = 0; i < sequence_length; ++i) {
#if 0
		std::cout << int(redueced_code_map[subsequence0[i]]) << ","
		<< int(redueced_code_map[subsequence1[i]]) << std::endl;
#endif
		mismatch_count +=
				(redueced_code_map[subsequence0[i]]
						!= redueced_code_map[subsequence1[i]]) ? 1 : 0;
	}
	return mismatch_count;
}

#endif /* SEED_SEARCHER_COMMON_H_ */
