# This algorithm finds all occurrences of a pattern in a text in linear time
# for preproccesing which help me then to construct bad char and good suffix matrices
# so i can say taht using it we can find l&r
# if():
# else:
# if():
# else:

"""
The bad character rule is a useful heuristic for mismatches near the right end of P,
but it has no effect
if the mismatching character from T occurs in P to the right of the mismatch point.
"""

# neither the good suffix rule nor the bad character rule shift P so far as
# to miss any occurrence of P. So the Boyer-Moore algorithm shifts by the largest amount
# given by either of the rules


# how much P should be "shifted? no shifting happening in the actual implementation.
# but the index k is increased to the point where the right end of P would be "shifted".
# Hence, each act of shifting P takes constant time.

# goes through from right to left
# and ofr know it is just like naive
# But GS &bchar are the diffrence makers


from collections import Counter


def zArray(s):
    n = len(s)
    try:
        n > 1
    except:
        print("impty string cant do it")

    zArr = [n] + [0] * (n - 1)

    right = 0
    left = 0
    for idx in range(1, n):
        try:
            zArr[idx] == 0
        except:
            print("error")
        # no prefix substring that starts before idx and ends after it

        if idx > right:
            n = 0
            while n + idx < len(s) and s[n] == s[n + idx]:
                n += 1
            zArr[idx] = n
            if n > 0:
                left = idx
                right = idx + n - 1
        else:
            p = idx - left

            if zArr[p] < right - idx + 1:
                zArr[idx] = zArr[p]
            else:
                i = right + 1
                while i < n and s[i] == s[i - idx]:
                    i += 1
                zArr[idx] = i - idx

                left = idx
                right = i - 1
    return zArr


# N array, maximum length of suffix of s
# which is also a suffix of s
# and cause it is for suffix i have to send to z array

def n_array(s):
    # it reverse s
    # the do preprocessing add it
    # reverse it again in the end
    return zArray(s[::-1])[::-1]


def leftPrimeArray(p, arr):
    # Compile array using p and  arr
    # L'[i] = largest index j less than n
    # such that N[j] = |P[i:]|
    n = len(p)
    lp = [0] * n
    idx = 0
    while (idx < (n - 1)):
        i = n - arr[idx]
        if (i < n):
            lp[i] = idx + 1
        idx = idx + 1
    return lp


def smallLeftPrimeArray(arr):
    n = len(arr)
    small_lp = [0] * n
    for i in range(n):
        if (arr[i] == i + 1):  # prefix matching a suffix
            small_lp[n - i - 1] = i + 1
    for i in range(n - 2, -1, -1):  # "smear" them out to the left
        if (small_lp[i] == 0):
            small_lp[i] = small_lp[i + 1]
    return small_lp


def goodSuffixTable(p):
    """ Return tables needed to apply good suffix rule. """
    n = n_array(p)
    lp = leftPrimeArray(p, n)
    # L[i] = largest index j less than n such that N[j] >= |P[i:]|
    l = [0] * len(p)
    l[1] = lp[1]
    i = 2
    while (i < len(p)):
        l[i] = max(l[i - 1], lp[i])
        i = i + 1

    return l, smallLeftPrimeArray(n)


def goodSuffixMismatch(i, bigLeftPrime, smallLeftPrime):
    """ Given a mismatch at offset i, and given L/L' and l' arrays,
        return amount to shift as determined by good suffix rule. """
    n = len(bigLeftPrime)
    if i == n - 1:
        return 0
    i = i + 1  # i points to leftmost matching position of P
    if bigLeftPrime[i] > 0:
        return n - bigLeftPrime[i]
    return n - smallLeftPrime[i]


# the right-to-left scan the bad character shift rule
# & the good suffix shift rule

def badCharTable(p, map):
    """ Given pattern string and list with ordered alphabet characters, create
        and return a dense bad character table.  Table is indexed by offset
        then by character. """
    tab = []
    nxt = [0] * len(map)
    i = 0
    while (i < len(p)):
        c = p[i]
        tab.append(nxt[:])
        nxt[map[c]] = i + 1
        i = i + 1
    return tab


class BoyerMoore():
    """ Encapsulates pattern and associated Boyer-Moore preprocessing. """

    # alphabet='ACGT' cause it is a DNA
    def __init__(self, p, alphabet='ACGT'):
        self.p = p
        self.alphabet = alphabet
        # Create map from alphabet characters to integers
        self.map = Counter(p)  # act as dict ofcourse it dict but much more effinect in freq problems {}
        for i in range(4):
            self.map[self.alphabet[i]] = i
        # Make bad character rule table
        self.bad_char = badCharTable(p, self.map)
        # Create good suffix rule table
        self.big_l, self.small_l_prime = goodSuffixTable(p)

    def bad_character_rule(self, i, c):
        """ Return # skips given by bad character rule at offset i """
        ci = self.map[c]
        return i - (self.bad_char[i][ci] - 1)

    def good_suffix_rule(self, i):
        # Given a mismatch at offset i, return amount to shift
        n = len(self.big_l)
        if i == n - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return n - self.big_l[i]
        return n - self.small_l_prime[i]

    def boyerMoore(self, p, t):
        """ Do Boyer-Moore matching """
        i = 0
        occurrences = []
        while i < len(t) - len(p) + 1:
            shift = 1
            mismatched = False
            for j in range(len(p) - 1, -1, -1):
                if p[j] != t[i + j]:
                    skip_bc = self.bad_character_rule(j, t[i + j])
                    skip_gs = self.good_suffix_rule(j)
                    shift = max(shift, skip_bc, skip_gs)
                    mismatched = True
                    break
            if not mismatched:
                occurrences.append(i)
                # special case of goodSuffix when p is exactly t
                skip_gs = len(self.small_l_prime) - self.small_l_prime[1]
                shift = max(shift, skip_gs)
            i += shift
        return occurrences


def testOne():
    t = 'GCTAGCTCTACGAGTCTA'
    p = 'CTT'
    p_bm = BoyerMoore(p, alphabet='ACGT')
    output = p_bm.boyerMoore(p, t)
    if not output:
        return -1
    else:
        return output


def testTwo():
    t = 'GCTAGCTCTACGAGTCTA'
    p = 'TCTA'
    p_bm = BoyerMoore(p, alphabet='ACGT')
    output = p_bm.boyerMoore(p, t)
    if not output:
        return -1
    else:
        return output


import unittest


class TestBoyerMoore(unittest.TestCase):
    def test_one(self):
        t = 'GCTAGCTCTACGAGTCTA'
        p = 'CTT'
        p_bm = BoyerMoore(p, alphabet='ACGT')
        output = p_bm.boyerMoore(p, t)
        self.assertFalse(output)

    def test_two(self):
        t = 'GCTAGCTCTACGAGTCTA'
        p = 'TCTA'
        p_bm = BoyerMoore(p, alphabet='ACGT')
        output = p_bm.boyerMoore(p, t)
        self.assertEqual(output, [6, 14])


def main():
    print(testOne())
    print(testTwo())


if __name__ == "__main__":
    main()
    unittest.main()
