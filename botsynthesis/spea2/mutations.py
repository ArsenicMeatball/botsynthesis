import copy
import logging
import random
import uuid

from Bio.Data import CodonTable

import botsynthesis.utils.dict_functions as dictf
from botsynthesis.utils.constants import MAX_POPULATION_SIZE, SEQUENCE_KEY


def mutate_codon(codon: str, codon_table: CodonTable):
    """Takes a codon and attempts to mutate it
    :param codon_table:
    :param codon:
    :return:
    """
    if len(codon) != 3:
        raise ValueError(
            "Codon {0} is not of size 3".format(codon)
        )  # not handled by biopython
    # get our amino acid to figure out what codon possibilities there are
    try:
        amino = codon_table.forward_table[codon]
    except KeyError as e:
        if codon in codon_table.stop_codons:
            logging.debug(
                f"Codon {codon} is a stop codon for codon table {codon_table}"
            )
            return codon
        logging.error(
            f"Codon {codon} is not handled in codon table {codon_table}"
        )
        raise e
    full_back_table = dictf.invert_dict(codon_table.forward_table)
    possible_codons = full_back_table[amino]
    if len(possible_codons) == 1:
        logging.info(
            "wasted a call on mutate codon - "
            "codon is unique and cannot be changed"
        )
        return codon
    # remove original codon
    possible_codons.remove(codon)
    logging.debug(
        "possible choices {}, initial codon {}, amino {}".format(
            possible_codons, codon, amino
        )
    )
    # pick another
    new_codon = random.choice(list(possible_codons))
    return new_codon


def mutate_seq(
        sequence: str,
        mutation_chance: float,
        codon_table: CodonTable,
) -> str:
    """
        mutate a DNA Seq based on the mutation probability,
        returns a different DNA Seq
    :param codon_table:
    :param mutation_chance:
    :param sequence:
    :return:
    """
    if len(sequence) < 3 or len(sequence) % 3 != 0:
        raise ValueError(
            "Sequence {} cannot be sliced into codons".format(sequence)
        )
    if mutation_chance < 0 or mutation_chance > 1:
        raise ValueError(
            "Mutation chance {} needs to be a float between 0 and 1".format(
                mutation_chance
            )
        )
    logging.debug(
        "Mutation % {0}, Codon Usage {1}, Unmutated Sequence {2}".format(
            mutation_chance, codon_table, sequence
        )
    )
    codons = []
    for idx in range(0, len(sequence), 3):
        # either add the og codon or mutated boy
        codon = str(sequence[idx: idx + 3])
        # random random generates random between 0 and 1,
        # if it is smaller than mutation chance, then mutate
        if random.random() < mutation_chance:
            codon = mutate_codon(codon, codon_table)
        codons.append(codon)
    # add everything at once because memory
    new_sequence = "".join(codons)
    logging.debug("Mutated sequence {0}".format(new_sequence))
    return new_sequence


def initialize_population(
        desired_population_size: int,
        parent_sequence: str,
        mutation_chance: float,
        codon_table: CodonTable,
) -> dict:
    """Initialize a population based around the parent sequence
    :param desired_population_size: the desired size of the population
    :param parent_sequence: the parent sequence - better being codon optimized
    :param mutation_chance: the chance to mutate each codon
    :param codon_table: dict mapping codons to amino acids. see CodonTable docs
    :return: dict containing the population
        {seq_id : {seq_key : actual sequence}, ...more seq ids}
    """
    if (
            desired_population_size < 1
            or MAX_POPULATION_SIZE < desired_population_size
    ):
        raise ValueError(
            "Desired population size {} is not within 1 and {}".format(
                desired_population_size, MAX_POPULATION_SIZE
            )
        )
    population = {}
    sequences = set()
    attempts = 0
    max_attempts = 25
    while (
        len(population) < desired_population_size and attempts < max_attempts
    ):
        seq = mutate_seq(parent_sequence, mutation_chance, codon_table)
        if seq not in sequences:
            sequences.add(seq)
            seq_id = uuid.uuid4()
            population[seq_id] = {SEQUENCE_KEY: seq}
        else:
            attempts += 1
            logging.debug('failed to create a "new" sequence')
    return population


def tournament_selection_without_replacement(
        population: dict,
        n_ary: int,
        minima: bool,
        fitness_key_name: str,
) -> str:
    """
    Get winner of a tournament (key of winner) based on the fitness value
    :param fitness_key_name:
    :param minima: look for lowest value if True(default),
    False looks for largest value
    :param n_ary: number of individuals to conduct tournament on,
    default = binary
    :param population: all the individuals we are selecting from
    :return: the key of the winner of tournament
    """
    if n_ary < 2:
        raise ValueError(
            f"tournament selection must be binary or larger, currently {n_ary}"
        )
    population_keys = list(population.keys())
    # choose a random one to start
    current_key = random.choice(population_keys)
    val_current = population[current_key][fitness_key_name]
    for _ in range(1, n_ary):
        contender_key = random.choice(population_keys)
        val_contender = population[contender_key][fitness_key_name]
        replace = (
            val_contender < val_current
            if minima
            else val_contender > val_current
        )
        if replace:
            current_key = contender_key
            val_current = val_contender
    # key that won tournament
    return current_key


def generate_mating_pool_from_archive(
        archive: dict,
        desired_mating_pool_size: int,
        fitness_key_name: str,
) -> dict:
    """Selects a portion of the archive to generate a mating pool
    :param fitness_key_name:
    :param archive:
    :param desired_mating_pool_size:
    :return:
    """
    # 1 < mating pool < archive
    if 1 > desired_mating_pool_size or len(archive) < desired_mating_pool_size:
        raise ValueError(
            f"Desired mating pool size {desired_mating_pool_size} "
            f"must be greater than 1, and must be smaller than "
            f"archive size {len(archive)}"
        )

    mating_pool = {}
    while len(mating_pool) < desired_mating_pool_size:
        # use a minima binary tournament based on the fitness
        # of the individual to add to mating pool
        key_to_add = tournament_selection_without_replacement(
            archive, n_ary=2, minima=True, fitness_key_name=fitness_key_name
        )
        mating_pool[key_to_add] = archive[key_to_add]
    return mating_pool


def recombine_dna_sequence(
        sequence1: str,
        sequence2: str,
        number_of_sites: int,
        spread_between_sites: int,
) -> str:
    """Recombines 2 sequences to try and create a new sequence incorporating
     portions of each
    :param sequence1:
    :param sequence2:
    :param number_of_sites:
    :param spread_between_sites:
    :return:
    """
    # recombination we must ensure at codon lengths only
    # (spread between sites = 3)
    if number_of_sites < 1:
        raise ValueError("number of sites needs to be greater than 0")
    if len(sequence1) != len(sequence2):
        raise ValueError("Cannot recombine sequences of different sizes")
    if sequence1 == sequence2:
        logging.warning("Cannot recombine the same sequence")
        return sequence1
    # randomly select all the crossover sites
    # get all possible sites and start chopping them off
    all_possible_crossover_sites = [
        num for num in range(0, len(sequence1), spread_between_sites)
    ]
    # pare down the crossover sites (may result in less than number of
    # desired sites but that is ok. RNGesus has spoken
    crossover_sites = {
        random.choice(all_possible_crossover_sites)
        for _ in range(number_of_sites)
    }
    crossover_sites = list(crossover_sites)
    crossover_sites.sort()
    # assuming sequences are randomly chosen, it is ok to pick 1 to go first
    current_sequence = sequence1
    last_pos = 0
    new_seq = ""
    chunks = []
    for site in crossover_sites:
        chunks.append(current_sequence[last_pos:site])
        current_sequence = (
            sequence1 if current_sequence == sequence2 else sequence2
        )
        last_pos = site
    # add the last bit
    chunks.append(current_sequence[last_pos:])
    # add them all together
    new_seq = new_seq.join(chunks)
    if len(new_seq) != len(sequence1):
        raise ValueError(
            f"Recombined sequence {len(new_seq)} is not the same size as "
            f"its ancestor {len(sequence1)}"
            f"\nBuilt from chunks {chunks}\n"
            f"Parents:\n{sequence1}\n{sequence2}"
        )
    return new_seq


def get_rec_sites_for_len(
        sequence_length: int, recombination_chance: float
) -> int:
    if recombination_chance > 1 or recombination_chance < 0:
        raise ValueError(
            "Recombination chance {} needs to be between 0 and 1".format(
                recombination_chance
            )
        )
    num_recomb_sites = round(sequence_length / 3 * recombination_chance)
    return num_recomb_sites


def generate_population_from_archive(
        _archive: dict,
        desired_mating_pool_size: int,
        num_recombination_sites: int,
        mutation_chance: float,
        desired_population_size: int,
        fitness_key_name: str,
        codon_table: CodonTable,
) -> dict:
    if desired_population_size < 1:
        raise ValueError(
            "desired population size ({}) must be greater than 0".format(
                desired_population_size
            )
        )
    # don't touch archive
    archive_copy = copy.deepcopy(_archive)
    # save all ids and sequences for quick look up
    archive_sequences_and_ids = {}
    for ide, attr in _archive.items():
        archive_sequences_and_ids[attr[SEQUENCE_KEY]] = ide

    # create mating pool as subset of copied archive
    mating_pool = generate_mating_pool_from_archive(
        archive_copy, desired_mating_pool_size, fitness_key_name
    )
    # try and create a population from the mating pool
    population = {}
    max_attempts = desired_population_size / 2
    # save newly generated sequences for quick lookup,
    # duplication will result in attempts
    sequences = set()
    attempts = 0
    while len(population) < desired_population_size:
        if attempts > max_attempts:
            logging.warning(
                f"Failed to create a complete "
                f"population of size {desired_population_size}, "
                f"only has {len(population)}"
            )
            break
        # get 2 sequences, attempt to recombine them into one, mutate it,
        # then add it in
        tourney_winner_key_1 = tournament_selection_without_replacement(
            population=mating_pool,
            n_ary=2,
            minima=True,
            fitness_key_name=fitness_key_name,
        )
        tourney_winner_key_2 = tournament_selection_without_replacement(
            population=mating_pool,
            n_ary=2,
            minima=True,
            fitness_key_name=fitness_key_name,
        )
        recombined_child = recombine_dna_sequence(
            sequence1=mating_pool[tourney_winner_key_1][SEQUENCE_KEY],
            sequence2=mating_pool[tourney_winner_key_2][SEQUENCE_KEY],
            number_of_sites=num_recombination_sites,
            spread_between_sites=3,
        )
        mutant_child = mutate_seq(
            sequence=recombined_child,
            mutation_chance=mutation_chance,
            codon_table=codon_table,
        )
        # check if we created a new one
        if mutant_child in archive_sequences_and_ids:
            # if it already exists then just add the og which
            # has a lot of information
            population[archive_sequences_and_ids[mutant_child]] = archive_copy[
                archive_sequences_and_ids[mutant_child]
            ]
            attempts += 1
        elif mutant_child in sequences:
            # failed to create a new sequence
            attempts += 1
        else:
            population[uuid.uuid4()] = {SEQUENCE_KEY: mutant_child}
            sequences.add(mutant_child)
    return population
