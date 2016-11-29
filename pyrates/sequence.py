"""Classes and functions to handle sequence data.
"""
import pyrates.utils as utils

class SequenceWithQuality(object):
    """A sequence and its quality scores.

    Args:
        sequence (:obj:`str`): The sequence to be stored.
        quality (:obj:`str`): Quality scores for base calls in sequence
            consisting of ASCII encoded phred scores.
        name (:obj:`str`, optional): Name to use for sequence.
    """
    __slots__ = '_sequence', '_quality', '_len', 'name'
    def __init__(self, sequence, quality, name=''):
        if len(sequence) != len(quality):
            raise ValueError("Sequence and quality have to have same length.")
        self._sequence = sequence
        self._quality = quality
        self._len = len(sequence)
        self.name = name

    @property
    def sequence(self):
        """Access to stored sequence :obj:`str`.

        Attempts to assign a new sequence that differs in length from
        the existing sequence will raise a :obj:`ValueError`.
        """
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        if len(value) != len(self):
            raise ValueError("Sequence of length %d expected, got %r (length: %d)." % \
                             (len(self), value, len(value)))
        self._sequence = value

    @property
    def quality(self):
        """Access to stored quality :obj:`str`.

        Attempts to assign a new qualities that differs in length from
        the existing ones will raise a :obj:`ValueError`.
        """
        return self._quality

    @quality.setter
    def quality(self, value):
        if len(value) != len(self):
            raise ValueError("Qualities of length %d expected, got %r (length: %d)" % \
                             (len(self), value, len(value)))
        self._quality = value

    def grosslydifferent(self, other, length=10, tolerance=7):
        """Partial comparison of sequences to determine whether
        they are substantially different.

        Args:
            other (:obj:`pyrates.sequence.SequenceWithQuality`): Sequence to compare to.
            length (:obj:`int`, optional): length of the prefix to test.
            tolerance (:obj:`int`, optional): Maximum number of mismatches allowed.

        Returns:
            :obj:`bool`: True if the sequences are considered too different
            to warrant further comparison.
        """
        diff = 0
        for i in range(length):
            if self._sequence[i] != other.sequence[i]:
                diff = diff + 1
        if diff > tolerance:
            return True
        return False

    def __str__(self):
        """Convert sequence + quality to FASTQ format."""
        return "@%s\n%s\n+\n%s" % (self.name, self.sequence, self.quality)

    def __repr__(self):
        return "SequenceWithQuality(sequence=%r, quality=%r, name=%r)" % \
                (self.sequence, self.quality, self.name)

    def __len__(self):
        return self._len

    def __getitem__(self, key):
        return (self._sequence[key], self._quality[key])

class SequenceStore(object):
    """Store a collection of sequences.

    This behaves like a set in the sense that each unique sequence is only represented
    once. Supports fast lookup of approximate matches.

    Args:
        max_length (:obj:`int`): Maximum sequence length supported by this store.
        alphabet (:obj:`tuple`): A list of all valid sequence characters.
        wildcard (:obj:`string`): A character that should be treated as a wildcard,
            i.e. match all letters in the alphabet.
    """
    __slots__ = '_alphabet', '_composition', '_index'
    _logger = utils.get_logger(__name__)

    def __init__(self, max_length, alphabet=('A', 'C', 'G', 'T')):
        self._index = {}
        self._alphabet = alphabet
        self._composition = {letter:[set() for _ in range(max_length+1)] for letter in alphabet}

    @classmethod
    def from_list(cls, sequences, **kw):
        """Create Sequence store from a list of sequences.

        Args:
            sequences (:obj:`list`): A list of sequences.

            Additional named arguments will be passed to the SequenceStore constructor.

        Returns:
            :obj:`pyrates.sequence.SequenceStore`
        """
        store = cls(len(sequences[0]), **kw)
        for seq in sequences:
            store.add(seq)
        return store

    def add(self, sequence, wildcard=None):
        """Add a sequence to the store.

        Args:
            sequence (:obj:`string`): New sequence to be added.
        """
        if sequence not in self._index:
            self._index[sequence] = {}
            if wildcard is None:
                for letter in self._alphabet:
                    letter_count = sequence.count(letter)
                    self._index[sequence][letter] = (letter_count, letter_count + 1)
                    self._composition[letter][letter_count].add(sequence)
            else:
                wilds = sequence.count(wildcard)
                for letter in self._alphabet:
                    letter_count = sequence.count(letter)
                    letter_index = (letter_count, letter_count + wilds + 1)
                    self._index[sequence][letter] = letter_index
                    for i in range(*letter_index):
                        self._composition[letter][i].add(sequence)

    def remove(self, item):
        """Remove a sequence from the sequence store.

        Args:
            item (:obj:`string`): Sequence to be removed.

        Raises:
            KeyError: if the sequence doesn't exist in the store.
        """
        for letter in self._alphabet:
            for i in range(*self._index[item][letter]):
                self._composition[letter][i].remove(item)
        del self._index[item]

    def discard(self, item):
        """Remove a sequence from the store if it exists.

        Args:
            item (:obj:`string`): Sequence to be removed.
        """
        if item in self._index:
            self.remove(item)

    def find(self, sequence, max_diff, wildcard=None):
        """Find best match for sequence in the store.

        Args:
            sequence (:obj:`string`): Sequence to search for.
            max_diff (:obj:`int`): Maximum number of mismatches allowed for a match.
            wildcard (:obj:`str`, optional): A character that should be treated as a wildcard.

        Returns:
            :obj:`tuple`: A tuple consisting of the best match found in the store and
            the number of differences between the returned match and the search string.
            If no suitable match was found `None` is returned instead.
        """
        match = self.search(sequence, 1, wildcard=wildcard)
        if len(match) == 0 or match[0][1] > max_diff:
            return None
        return match[0]

    def search(self, sequence, max_diff, max_hits=10, raw=False, wildcard=None):
        """Search the sequence store for all approximate matches to a search pattern.

        Args:
            sequence (:obj:`string`): Sequence to search for.
            max_hits (:obj:`int`, optional): Maximum number of results to return.
                set to _None_ to return all candidates. Ignored if `raw` is _True_.
            raw (:obj:`bool`, optional): Flag indicating whether the raw sequence
                matches should be returned instead of sequence/distance pairs.
            wildcard (:obj:`str`, optional): A character that should be treated as a wildcard.

        Returns:
            If `raw` is _True_ an unordered :obj:`list` of candidates is returned,
            otherwise a list of (sequence, distance) tuples is returned.
        """
        if sequence in self._index:
            if raw:
                return [sequence]
            return [(sequence, 0)]
        wilds = 0
        if wildcard is not None:
            wilds = sequence.count(wildcard)
        letter = self._alphabet[0]
        letter_count = sequence.count(letter)
        min_count = max(0, letter_count - max_diff)
        max_count = min(len(self._composition[letter]), letter_count + max_diff + wilds + 1)
        candidates = []
        for cand in self._composition[letter][min_count:max_count]:
            candidates.extend(cand)
        for letter in self._alphabet[1:]:
            letter_count = sequence.count(letter)
            min_count = max(0, letter_count - max_diff)
            max_count = min(len(self._composition[letter]), letter_count + max_diff + wilds + 1)
            for cand in self._composition[letter][min_count:max_count]:
                candidates.extend(cand)
        candidates = set(candidates)
        if raw:
            return candidates
        candidates = [(cand, self.diff(sequence, cand)) for cand in candidates]
        candidates.sort(key=lambda x: x[1])
        candidates = [cand for cand in candidates if cand[1] <= max_diff]
        if max_hits is not None:
            candidates = candidates[:max_hits]
        return candidates

    @staticmethod
    def diff(seq1, seq2):
        """Compute Hamming distance between two sequences.
        """
        diff = 0
        for (i, letter) in enumerate(seq1):
            if letter != seq2[i]:
                diff += 1
        return diff

    def __len__(self):
        return len(self._index)

    def __contains__(self, item):
        return item in self._index
