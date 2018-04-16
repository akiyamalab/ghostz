LAST_CC=g++
CC=gcc
CXX=g++
COMMON_FLAGS=-Wall -Wextra -pedantic -fopenmp -O3 -g
CC_FLAGS=${COMMON_FLAGS} 
LDLIBS	=-lm -fopenmp

ifeq ($(PROFILE), Yes)
CC_FLAGS=${COMMON_FLAGS} -g
endif


C_SRC=./ext/karlin/src/karlin.c 

CPP_SRC=./ext/seg/src/seg.cpp \
	./src/align_main.cpp \
	./src/aligner.cpp \
	./src/alphabet_coder.cpp \
	./src/alphabet_type.cpp \
	./src/chain_filter.cpp \
	./src/database.cpp \
	./src/database_build_main.cpp \
	./src/database_chunk.cpp \
	./src/dna_sequence.cpp \
	./src/dna_type.cpp \
	./src/edit_blocks.cpp \
	./src/fasta_sequence_reader.cpp \
	./src/gapped_extender.cpp \
	./src/k_mer_sequences_index.cpp \
	./src/main.cpp \
	./src/one_mismatch_hash_generator.cpp \
	./src/protein_query.cpp \
	./src/protein_sequence.cpp \
	./src/protein_type.cpp \
	./src/queries.cpp \
	./src/query.cpp \
	./src/reduced_alphabet_coder.cpp \
	./src/reduced_alphabet_file_reader.cpp \
	./src/reduced_alphabet_k_mer_hash_function.cpp \
	./src/reduced_alphabet_variable_hash_function.cpp \
	./src/score_matrix.cpp \
	./src/score_matrix_reader.cpp \
	./src/seed_searcher.cpp \
	./src/seed_searcher_common.cpp \
	./src/seed_searcher_database_parameters.cpp \
	./src/seed_searcher_query_parameters.cpp \
	./src/sequence.cpp \
	./src/sequence_no_filter.cpp \
	./src/sequence_seg_filter.cpp \
	./src/statistics.cpp \
	./src/translated_dna_query.cpp \
	./src/translator.cpp \
	./src/ungapped_extender.cpp \
	./src/variable_hash_clustering_seuences_index.cpp \
	./src/variable_hash_sequences_index.cpp 


OBJS =
OBJS += $(C_SRC:%.c=%.o)
OBJS += $(CPP_SRC:%.cpp=%.o)

.SUFFIXES:	.o

.PHONY: all
all:ghostz

ghostz: $(OBJS)
	$(LAST_CC) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)

.c.o:
	$(CC) -c $(CC_FLAGS) $< -o $@  $(INCLUDES)

.cpp.o:
	$(CXX) -c $(CC_FLAGS) $< -o $@  $(INCLUDES)



.PHONY: clean
clean:
	rm -f $(OBJS)
