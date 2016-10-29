"""Classes and functions to handle sequence data.
"""

class SequenceWithQuality(object):
    """A sequence and its quality scores.

    Args:
        sequence (:obj:`str`): The sequence to be stored.
        quality (:obj:`str`): Quality scores for base calls in sequence
            consisting of ASCII encoded phred scores.
        name (:obj:`str`, optional): Name to use for sequence.
    """
    __slots__ = '_sequence', '_quality', 'name'
    def __init__(self, sequence, quality, name=''):
        if len(sequence) != len(quality):
            raise ValueError("Sequence and quality have to have same length.")
        self._sequence = sequence
        self._quality = quality
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
        return len(self.sequence)

    def __getitem__(self, key):
        return (self._sequence[key], self._quality[key])
