"""Test sequence module."""

from nose2.tools import params
from nose2.tools.such import helper
from pyrates.sequence import SequenceWithQuality, SequenceStore, GroupedSequenceStore

def test_swq_new():
    """Create sequence objects"""
    read = "ACTGGGTGTAT"
    qual = "IIIIIIIIIII"
    name = 'test'
    seq = SequenceWithQuality(read, qual)
    seq = SequenceWithQuality(read, qual, name)
    assert len(seq) == len(read), "Lengths of sequence and read don't match"
    assert seq.sequence == read, "Sequence doesn't match input read"
    assert seq.quality == qual, "Quality doesn't match input"
    assert seq.name == name, "Name doesn't match input"

def test_swq_access():
    """Access to sequence object members"""
    read = "ACTGGGTGTAT"
    read2 = "TATGTGGGTCA"
    qual = "IIIIIIIIIII"
    qual2 = "DDDDDDDDDDD"
    name = 'test'
    seq = SequenceWithQuality(read, qual, name)
    assert seq[0] == ('A', 'I'), "Failed to extract sequence/quality pair"
    seq.sequence = read2
    assert seq.sequence == read2, "Failed to change read sequence"
    seq.quality = qual2
    assert seq.quality == qual2, "Failed to change qualities"

    ## check that attempts to set members to new values that
    ## would result in an inconsistent state raise an exception
    with helper.assertRaises(ValueError):
        SequenceWithQuality("ACTG", qual)
    with helper.assertRaises(ValueError):
        seq.sequence = "ACTG"
    with helper.assertRaises(ValueError):
        SequenceWithQuality(read, "!!!!")
    with helper.assertRaises(ValueError):
        seq.quality = "!!!!!"

def test_swq_str():
    """String representation of sequence objects"""
    read = "ACTGGGTGTAT"
    qual = "IIIIIIIIIII"
    name = 'test'
    seq = SequenceWithQuality(read, qual, name)
    expect_str = "@test\nACTGGGTGTAT\n+\nIIIIIIIIIII"
    expect_repr = "SequenceWithQuality(sequence='ACTGGGTGTAT', quality='IIIIIIIIIII', name='test')"
    assert str(seq) == expect_str, "%r != %r" % (str(seq), expect_str)
    assert repr(seq) == expect_repr, "%r != %r" % (repr(seq), expect_repr)

def test_grosslydifferent_equal():
    """Recognise identical sequences as similar"""
    seq1 = SequenceWithQuality('AACTGTGAGTGTAGATGTTC', 'I'*20)
    assert not seq1.grosslydifferent(seq1), \
        "Sequence %r considered too different from itself" % seq1.sequence

@params(('AACTGTGAGTGTAGATGTTC', 'AACTTTGAGTGTAGATGTTC'), \
        ('AACTGTGAGTGTAGATGTTC', 'TTTTTTTTTTGTAGATGTTC'))
def test_grosslydifferent_similar(seq1, seq2):
    """Tolerate small differences between sequences."""
    seq1 = SequenceWithQuality(seq1, 'I'*len(seq1))
    seq2 = SequenceWithQuality(seq2, 'I'*len(seq2))
    assert not seq1.grosslydifferent(seq2), \
        "Sequences %r and %r considered too different" % (seq1.sequence, seq2.sequence)

def test_grosslydifferent_diff():
    """Tolerate small differences between sequences."""
    seq1 = SequenceWithQuality('AACTGTGAGTGTAGATGTTC', 'I'*20)
    seq2 = SequenceWithQuality('GTAGATGTTCGTAGATGTTC', 'I'*20)
    assert seq1.grosslydifferent(seq2), \
        "Sequences %s and %s considered similar" % (seq1.sequence, seq2.sequence)

def test_store_new():
    """Create empty sequence store"""
    store = SequenceStore(5)
    assert len(store) == 0, "%r != 0" % len(store)

def test_store_add():
    """Add entries to sequence store"""
    store = SequenceStore(4)
    store.add("AAAA")
    assert len(store) == 1, "%r != 1" % len(store)
    assert "AAAA" in store
    store.add("CTGT")
    assert len(store) == 2, "%r != 2" % len(store)
    assert "CTGT" in store

def test_store_add_wild():
    """Add entries to sequence store"""
    store = SequenceStore(4)
    store.add("AAAA", wildcard='N')
    assert len(store) == 1, "%r != 1" % len(store)
    assert "AAAA" in store
    store.add("CTNT", wildcard='N')
    assert len(store) == 2, "%r != 2" % len(store)
    assert "CTNT" in store

def test_store_contains():
    """Test sequences for membership in store"""
    store = SequenceStore(4)
    store.add("AAAA")
    assert "AAAA" in store, "'AAAA' not found in store"
    assert "AACT" not in store, "'AACT' should not be in store"

def test_store_remove():
    """Remove sequences from store"""
    store = SequenceStore(4)
    store.add("AAAA")
    store.add("CTGT")
    assert "AAAA" in store, "'AAAA' not found in store"
    assert "CTGT" in store, "'CTGT' not found in store"
    store.remove("AAAA")
    assert "AAAA" not in store, "'AAAA' remains in store after removal"
    with helper.assertRaises(KeyError):
        store.remove("AAAA")
    assert "CTGT" in store, "'CTGT' not found in store"
    assert len(store) == 1, "%r != 1" % len(store)
    store.remove("CTGT")
    assert "CTGT" not in store, "'CTGT' remains in store after removal"
    assert len(store) == 0, "%r != 0" % len(store)

def test_store_discard():
    """Discard sequences from store if they exist"""
    store = SequenceStore(4)
    store.add("AAAA")
    store.add("CTGT")
    assert "AAAA" in store, "'AAAA' not found in store"
    assert "CTGT" in store, "'CTGT' not found in store"
    store.discard("AAAA")
    assert "AAAA" not in store, "'AAAA' remains in store after removal"
    store.discard("AAAA")
    assert "CTGT" in store, "'CTGT' not found in store"
    assert len(store) == 1, "%r != 1" % len(store)
    store.discard("CTGT")
    assert "CTGT" not in store, "'CTGT' remains in store after removal"
    assert len(store) == 0, "%r != 0" % len(store)

def test_store_search():
    """Find all approximate matches"""
    store = SequenceStore(4)
    store.add("AAAA")
    store.add("AAAT")
    store.add("AATT")
    store.add("ATTT")
    match = store.search('TTTT', 4, max_hits=None)
    assert len(match) == 4, "%r != 4" % len(match)
    match = store.search('TTTT', 2, max_hits=None)
    assert len(match) == 2, "%r != 2" % len(match)


@params(('AAAA', ('AAAA', 0)), ('CATT', ('AATT', 1)), ('GGGG', None))
def test_store_find(search, expect):
    """Find best approximate match"""
    store = SequenceStore(4)
    store.add("AAAA")
    store.add("AATT")
    store.add("TTTT")
    match = store.find(search, 2)
    assert match == expect, "%r != %r" % (match, expect)

@params(('AANA', ('AAAA', 1)), ('CATT', ('AATT', 1)), ('NGGG', None))
def test_store_find_wild(search, expect):
    """Find best approximate match"""
    store = SequenceStore(4)
    store.add("AAAA")
    store.add("AATT")
    store.add("TTTT")
    match = store.find(search, 2, wildcard='N')
    assert match == expect, "%r != %r" % (match, expect)

def test_store_list():
    """Create sequence store from list."""
    store = SequenceStore.from_list(['AAAA', 'CTGT'])
    assert "AAAA" in store, "'AAAA' not found in store"
    assert "CTGT" in store, "'CTGT' not found in store"

def test_grouped_new():
    """Create GroupedSequenceStore"""
    store = GroupedSequenceStore(4, tag_size=2)
    assert len(store) == 0, "%r != 0" % len(store)
    assert "AAAA" not in store

def test_grouped_add():
    """Add entries to GroupedSequenceStore"""
    store = GroupedSequenceStore(4, tag_size=2)
    store.add("AAAA")
    assert len(store) == 1, "%r != 1" % len(store)
    assert "AAAA" in store
    store.add("CTGT")
    assert len(store) == 2, "%r != 2" % len(store)
    assert "CTGT" in store
    store.add("CTAA")
    assert len(store) == 3, "%r != 3" % len(store)
    assert "CTAA" in store

def test_grouped_add_wild():
    """Add entries to GroupedSequenceStore"""
    store = GroupedSequenceStore(4, tag_size=2, wildcard='N')
    store.add("AAAA")
    assert len(store) == 1, "%r != 1" % len(store)
    assert "AAAA" in store
    store.add("CTNT")
    assert len(store) == 2, "%r != 2" % len(store)
    assert "CTNT" in store, "'CTNT' not found in store"

def test_grouped_contains():
    """Test sequences for membership in GroupedSequenceStore"""
    store = GroupedSequenceStore(4, tag_size=2, wildcard='N')
    store.add("AAAA")
    store.add("NAAA")
    assert "AAAA" in store, "'AAAA' not found in store"
    assert "NAAA" in store, "'NAAA' not found in store"
    assert "AACT" not in store, "'AACT' should not be in store"
    assert "AANT" not in store, "'AANT' should not be in store"

def test_grouped_remove():
    """Remove sequences from GroupedSequenceStore"""
    store = GroupedSequenceStore(4, tag_size=2)
    store.add("AAAA")
    store.add("CTGT")
    assert "AAAA" in store, "'AAAA' not found in store"
    assert "CTGT" in store, "'CTGT' not found in store"
    store.remove("AAAA")
    assert "AAAA" not in store, "'AAAA' remains in store after removal"
    with helper.assertRaises(KeyError):
        store.remove("AAAA")
    assert "CTGT" in store, "'CTGT' not found in store"
    assert len(store) == 1, "%r != 1" % len(store)
    store.remove("CTGT")
    assert "CTGT" not in store, "'CTGT' remains in store after removal"
    assert len(store) == 0, "%r != 0" % len(store)

def test_grouped_discard():
    """Discard sequences from GroupedSequenceStore if they exist"""
    store = GroupedSequenceStore(4, tag_size=2)
    store.add("AAAA")
    store.add("CTGT")
    assert "AAAA" in store, "'AAAA' not found in store"
    assert "CTGT" in store, "'CTGT' not found in store"
    store.discard("AAAA")
    assert "AAAA" not in store, "'AAAA' remains in store after removal"
    store.discard("AAAA")
    assert "CTGT" in store, "'CTGT' not found in store"
    assert len(store) == 1, "%r != 1" % len(store)
    store.discard("CTGT")
    assert "CTGT" not in store, "'CTGT' remains in store after removal"
    assert len(store) == 0, "%r != 0" % len(store)

def test_grouped_search():
    """Find all approximate matches"""
    store = GroupedSequenceStore(4, tag_size=2)
    store.add("AAAA")
    store.add("AAAT")
    store.add("AATT")
    store.add("ATTT")
    match = store.search('TTTT', 4, max_hits=None)
    assert len(match) == 2, "%r != 2 (found %r)" % (len(match), match)
    match = store.search('TTTT', 2, max_hits=None)
    assert len(match) == 2, "%r != 2 (found %r)" % (len(match), match)

@params(('AAAA', ('AAAA', 0)), ('CATT', ('AATT', 1)), ('GGGG', None))
def test_grouped_find(search, expect):
    """Find best approximate match"""
    store = GroupedSequenceStore(4, tag_size=2)
    store.add("AAAA")
    store.add("AATT")
    store.add("TTTT")
    match = store.find(search, 2)
    assert match == expect, "%r != %r" % (match, expect)

@params(('AANA', ('AAAA', 1)), ('CATT', ('AATT', 1)), ('NGGG', None))
def test_grouped_find_wild(search, expect):
    """Find best approximate match"""
    store = GroupedSequenceStore(4, tag_size=2, wildcard='N')
    store.add("AAAA")
    store.add("AATT")
    store.add("TTTT")
    match = store.find(search, 2)
    assert match == expect, "%r != %r" % (match, expect)
