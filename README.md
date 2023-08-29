# bsplit-py
Bernoulli Split (python)

# Summary

This python library is a proving ground for calculating the necessary coefficients for a novel algorithm called BSplit 
(Bernoulli Split).

This library is intended as a pre-cursor to the bsplit-rs (rust) project.

This library can produce the various coefficients and rust-syntax for arrays that would be tedeous to do by hand.

# Bernoulli

The bernoulli Distribution plots the ratio of the number of possible combinations that a binary system can have 
true and false pairs. When almost all are True OR when almost all are False, then there are very few possible 
combinations.  When there is a relatively even number of True and False, the quantity grows semi-factorally.  
This produces something LIKE a bell-shaped curve.

[Bernoulli Distribution](https://en.wikipedia.org/wiki/Bernoulli_distribution)

![Bernoulli Graph](https://upload.wikimedia.org/wikipedia/commons/3/3a/Standard_deviation_diagram_micro.svg)

This project contains my various explorations over the years of finding mathmatical boundaries for use in 
search/discovery and compression (with random access)

# Sterling/ Shannon

The bernoulli distribution is the basis of the Shannon Entropy approximation:

```python
def shannon_entropy(probs:list[float])->float: 
    return -sum([ p * math.log2(1/p) for p in probs]) 
```

$$ShannonEntropy = - \sum_{i=0}^{K} p_i * \log_2 p_i  $$

This equation is only an approximation because:
For a binary system p and 1.0-p the probability is irrational, yet the actual ratio is rational (two whole numbers).  
Namely if I have 32-choose-5, there is a discrete number of possible combinations, yet the shannon equation would yield:

```python
x_estimaged=shannon_entropy([5/32, (32-5)/32]) # 0.6
yet the actual combinations:
def entropy(n:int, k:int)->float:
    return math.log2(math.comb(n,k))
x_actual=entropy(32,5) / 32 # 0.55061
```

In other words, treating a binary system of `p = true_samples / total_samples` is wildly innacurate until you start to 
approach `total_samples = infinity`.

For VERY small `total_samples`, quantization is so dominant that the error can be thousands of percent off.

Yet, the sterling-approximation, which is used as the basis for the Shannon-approximation has been the gold standard 
for all topics concerning entropy.

# BSplit Estimator

An alternate derivation of the sterling-approximation provides a simpler binary equation:

```python
def clg(n:int)->int:
    return math.ceil(math.log2(n))
def bsplit_size(n:int,k:int)->int:
    return math.ceil( k * (clg(n)-clg(k+1) + 1.5) )
```

or rather

**bsplit-entropy calcuation**
$$entropy = k * (\log_2 n - \log_2 (k+1) + 1.5)$$

This is not a HARD limit (like Shannon), but provides a very very close upper-limit to the TRUE entropy level 
(defined as math.comb(n,k)), and through a specific algorithm can be shown to always only need 1.5 bits extra PER element (k).

This estimate is also not accurate within `1/4N < K < 3/4N`.

# ALGORITHM

* Separate a key-space into two parts; high-frequency and low-frequency
* Define a log-base-2 for the following

```python
    bn=clg(n)   # Nominal Bits per element
    bk=clg(k+1) # A pallet proportionate to log(k)
    br=bn-bk    # Residual bits
```

* The high-frequency has a size of `k * br` bits (exactly).
* You can random-access skip over the high-frequency bits because it is deterministic.. EACH entry K can be found at offset `start + k * br`
* Now store a low-frequency *pallet* bitmap size exactly `2**bk`.  This contains 1 bit for each partition of 2**bn.. If 1, there is at least 1 element (k)
* Now take the population of the pallet

```python
    u=pallet.bit_count()  # Number of 1s; number of partitions that have at least 1 element; u standard for 'unique'
    c=k-u                 # Number of collisions; e.g. number of elements that must be bunched up in another pallet partition
    collision_max_size = clg(u)*c-lookup(u,c)          # how many bits needed to "find" collisions (method-B)
    extra_bits = min(u - 1, collsion_max_size) # lesser of method-A and method-B (see below)
```

* This is the end of the 'fixed-length' portion; The deterministic portion.. The rest depends on the distribution.  What follows is the collision resolution bit-sequence (finite state machine).  There are two sub-flavors
    * Method-A: if c <= extra_bits: then use exactly 1 bit per "true" pallet bit minus one. Note we can infer the LAST partition, as it's whatever is left over, thus we have u-1
    * Method-B: else: then store a FiniteStateMachine who's upper-bound is `clg(u)*c`, but can save some bits in certain permutations

* At this point, we can easily decode the bitmap as an iterable bit-sequence.  In the case of Method-A, we can also easily decode the collision-bit-sequence.  For Method-B, we need to perform a decode of the second sequence to turn into a larger bit-sequence (of size u).
* If the goal is to SKIP all the bits in this encoding, we have all the information we need (after calculating pallet size and perform a bit-count).  Simply adjust the cursor past the end of the entire tripplet
* If the goal is to decode ALL the items, then simply activate a pair of bit-sequence iterators (generators) and walk

```python
    collision_iter=iter(extra_bits_slice)
    residual_iter=iter(residual_bit_slice, br)
    for header_prefix in iter(pallet):
        header = header_prefix << br
        if collision_iter.next():
            yield decode_group(residual_iter, header)
        else:
            yield header | residual_iter.next()
```

For items with no collisions, we simply yield br bits (prefixed with the header).

For items with collisions, we need to recursively decode a bundle of bits.

BUT, we don't know how many items have collided; all we know is there is AT LEAST 2.

Thus we can use the sub algorithm:

```python
cloned_iter=residual_iter.clone()
word0 = cloned_iter.next()
count = 2 + unitary_scan(word0, max_length=br) # return number from 0 .. br; thus 2 .. br + 2
if count > br:
    count = int(cloned_iter.next()) + br  # use second word AS the count; this is guaranteed to be 0..2**br
    # skip 2 full words 
    residual_iter.skip(2)
    encoded_bits = residual_iter.slice(count-2)
else:
    encoded_bits = residual_iter.slice(count) # consume and pack into a new bitmap count words (size br each)
    encoded_bits.skip(count-1) # skip over the unitary bits we used above
decode(2**br, count, encoded_bits) # use this slice and yield their values
```

## Edge cases
It can be shown that the above sequence never exceeds the bit size for any count of elements. 
In fact, if the density of the collisions is high enough; this wastes bits.. 
A non-random access version of this algorithm could have recovered them.


# Higher orders

Having a bitmap of size 1 million isn't necessarily that efficient; 
especially since you need to make the population count. It also hurts random-access performance.  
Thus using a chunked recursion of a machine-word-size is most optimal

```python
def recurse(bn, k, bits):
    # edge cases needed to be handled here

    # main route decision tree
    if k > 64:
       residual, pallet, collisions = split(bits)
       ki = [x for x in recurse(k+64-1, 63, bits)]
       br = bn-clg(k)
       for k in ki:
         yield recurse(br, k, bits)
    else:
       yield leaf_algorithm(br, k, bits)
```
# Alternate methods

In addition to yielding each entry (e.g. a FULL decode), we can binary-search (more accurate batch-search).

## Find ith item
This can be used to represent an offset+length map. Each number stored is the "address" in a file of the ending byte.
Thus to see an offset+length, you would 


```python
offset=findIth(i-1)
length=findIth(i)-offset
```

* skip sub-groups while a cumulative `k_left+ki` is less than `k`
* exit once the entry is found
* in leaf, skip any sub-partition where `k_left+collsion_count` is less than `k`

## Return True if address A found
This is the typical bitmap representing a Set.  You want to test if a number is in the set

* skip sub-groups while `n_left+2**bn` is less than target A
* exit once the entry is passed (whether found or not)
* in leaf, if pallet containing A is 0, return False
* in leaf, decode the collision-group containing A.  Walk each colliding item and return True if the suffix (br bits) matches. 
* Else return False

## Return index of set A or set B
This would be if I had two categories A or B.  You add N items, but either to list A or list B.  You then maintain a
bitmap: 0 for list A, 1 for list B of where it was routed. Now you don't need to store an item-type; the bitmap is 
sufficient. It also encodes WHERE the item is stored into it's respective list.

```python
index, is_in_left = findIndex(item_number)
if is_in_left:
    return left[index]
else:
    return right[index]
```

We use the same algorithm as `Return True if address A found` but instead of JUST returning a bool, we return n_left

Note this algorithm can be repeated.. On first repetition, you have 4-classifications:
left_left, left_right, right_left, right_right

On second repetition, you have 8 classifications, and so on.  If the ratio of items is similar, then this does not provide
much advantage over just using a small enumerate-integer (say a 4 bit int in the case of 16 classifications).  However, if
the ratio is 2:1 or ideally 10:1, then this will follow Shannons' equation of entropy (above).

# Special-case groupings

The above algorithm works well for between 8 and 64 elements. 
For larger, the recursive method produces less than optimal bit-sizes but allows for random access AND 
allows for 64bit operations (being computationally more efficient). It also works well for having BYTE-ALIGNED 
groupings, (such that byte-aligned operations stay on the left of the byte-array and bit-aligned are pruned 
from the right).

For entropy ratios of  1/4 < k/n < 3/4, however, it is better to use a constrained bitmap.  
Namely store the N-41 bits, then use arithmatic encoding to store an integer of max size `math.comb(41,k)` in the 
suffix.  
This guarantees at least 3 bits are removed from the bitmap (as the largest integer for 41-choose-20 is only 38 bits.  
41 is the smallest number N where this is true.  
If N is less than 41, then store N-9 bits, then use arithmatic encoding to store an integer 
of max size `math.comb(9,k)` in the suffix.  This guarantees at least 2 bits are removed 
(as 9-choose-4 requires 7 bits).  
If N is less than 9, then simply store the N-1 bits, and infer the last bit from the population count of the 
first N-1 bits. (e.g. similar to a parity encoding)

For K sizes > 3/4N, store the INVERSE of the bitmap.. 

```python
if 2*k > 2**bn: return inverse(recurse(bn,k-2**bn,bits))
```

For K sizes 0 and 1, N-1 and N, simply store 0-bits or clg(n) bits.  No need for a complex multi-level decoding.

For K sizes 2 .. 8 and N-8 .. N-2, store a finite-state-machine which removes an upper-bound of factorial(k) bits.  
The general approach is to "Split" N into two halves.  Then use a unitary bitmap to say what the imbalance is.  
Smaller bit-sequences are closer to even left/right distribution.  
Larger bit-sequences are heavily skewed left or right (requiring an extra bit to say which side is larger).  
Depending on the exact code-book chosen, it is possible to guarnatee that:

* n-choose-2 saves 1 bit (floor(log2(factorial(2))) == 1)
* n-choose-3 saves 2 bits (floor(log2(factorial(3))) == 2)
* n-choose-4 saves 4 bits (floor(log2(factorial(3))) == 4)
* n-choose-5 saves 6 bits (floor(log2(factorial(3))) == 6)
* n-choose-6 saves 8 bits (while fact(6) would have saved 9 bits, I've not found a pure F.S.M. that accomplishes this without massive tables)
* n-choose-7 saves 11 bits
* n-choose-8 saves 13 bits

The above numbers can be used to calculate the "skip" length. So, for example, bn=20, k=7, the total bits is `20*7-11` or 129 bits exactly.

Note, in the F.S.M. some combinations can terminate in fewer bits.. To allow for random access, the bits should FIRST be sliced over (thus wasting some potential compression).

The nature of the FSM is to recurse to specialty functions, So for example, in n_choose_5, we'd have:

```python
def n_choose_5(n, word):
    left_side=word.next(bits=1)==1
    depth=word.unitary_scan(max=2)
    k_left=2-depth if left_side else 3+depth # note, can use zigzag encoding to avoid the if-statement
    k_right=5-k_left
    k_small=min(k_left, k_right)
    n_left=n//2
    n_right=n-n_left # handle order swap
    n_small, n_large = ...
    match k_small:
        case 2: n_choose_2(n_small,word) + n_choose_3(n_large,word)
        case 1: n_choose_1(n_small,word) + n_choose_4(n_large,word)
        case 0: n_choose_5(n_large,word)
    # handle re-ordering
    # emit
```        
A CUSTOM F.S.M. is necessary for each combination. Odds tend to be easier than Evens. But odds above 5 add to the complexity BECAUSE of the Bernoulli distribution and the factorial-growth rate as we have larger and larger distribution-group-sizes (moving towards the center of the bernoulli curve). Hense the use of the prior algorithm to simplify the need to create custom FSMs for a large number of bernoulli arrangements.

Note, it is possible to construct F.S.M.s for arbitrarily large combinations, however, those in the range 1/4N < K < 3/4N become unstable, where some combinations contract to only K bits, or explode to 2N bits.  Thus the use of a separate encoding algorithm in that range is necessary.  Note it is also possible to use range-encoding to guarantee math.ceil(math.log2(math.comb(n,k))) bits.  Indeed this is used in 41-choose-K and 9-choose-K.  The number of tables necessary becomes prohibatively expensive, however.  Those tables do not reduce the bit-size by much v.s. a n-choose-k outside the 1/4th .. 3/4th range.

Note A GOAL is to allow **GPU-optimized** execution, thus the avoidance of if-statements where possible.

# Demonstration of entropy upper-bound

* An upper bound of `1 / 4 n  < k < 3 / 4 n`

Empirically, n-choose-n/2 (worst case, peak of bernoulli chart) is `n-1/2 log_2(n)` bits.  E.g. at n=1024, then 1024-choose-512
is ROUGHLY 1024-5 bits (or 1019). The actual value is 1018.67 bits, rounded up is 1019.

Thus at 128, we'd expect to save 3 bits; and that's achieved in our above algorithm that removes last 41 bits.

Starting at 256, we SHOULD be able to save 4 bits v.s. a raw bitmap, but, the computational complexity of arithmatically
enumerating 256-choose-128 is probably not worth saving a single bit out of 125 bits. (less than 1% savings)


* A trend for smaller K

For  `16 <= k <= 64` and `n==2**16` we have:

```python
pallet=64             # 64 bits is 2**clk(64)
max_collision_size=32 # not investigated here, edge cases can exceed this slightly
br=bn-clg(64)         # 16-6  => 10 bits per residual element
residue_size=k*br     # 64*10 => 640 bits
sum=64+32+640         # 736 bits
avg=sum/64            # 11.5 bits per element
est=k*(bn-clg(k+1)+1.5) # about to remove the k to produce an average
avg_est=est/k         # bn-clg(k+1)+1.5 => 16-6+1.5 => 11.5 bits per element
```

So notice our calculated estimate is roughly identical to the decomposed average.

We note that we make the assumption that max_collision_size is only half of K. This is contingent on being careful
with edge cases, and for performance reasons we don't do so. Instead, we ALLOW a variation of max_collision_size and require
our skip-over logic to decode the pallet and see the exact max_collision_size for each sub-group that needs skipping over.

* A trend for N-choose-8 and below

The FiniteStateMachine can be fully fleshed out for every possible permutation of N as a power of 2.

This can be accomplished by working out all Ns that are 1/4th N to 3/4th N which trigger most of the edge cases.
Thus we can work out any `N<=64` with all permutations of `K=2..8`.  Once they are verified, all sizes `N > 64` and still
a power of 2 will linearly add bits.

This process reveals that N-choose-6 SHOULD calculate a saving of 9.577 bits per element, but we are only able to
extract 8 bits per element in the worst case.  This is because the left/right split is unstable - being OVERLY compressive
for highly imbalanced distributions. This needs to average out; and that happens towards the uniformly distributed cases.

For example take a finite-state-machine for 1024-choose-6

```
k_left / k_right : header_bits left_bits right_bits => total_bits
0/6: 4 0  46 => 50
1/5: 4 9  39 => 52
2/4: 2 17 32 => 51
3/3: 2 25 25 => 52
4/2: 2 32 17 => 51
5/1: 4 39 9  => 52
6/0: 4 46 0  => 50
```

This compared to the actual entropy of `50.487` bits (rounded to 51). We see that our 0/6 and 6/0 left/right
distributions do BETTER than true-entropy, thus something somewhere has to do slightly worse. In fact, 0/6 can do
MUCH better than 50 as this is recursive. That assymetric amplifies both the best case (to be better) and worst case
(to be worse).


# Use cases A.I.

Nvidia utilizes a sparse-float method in **cuBLAS**.  Their method is to take a group of 4 coefficients
and select exactly 2. The other 2 are treated as zero.
This allows
1. Halving the number of multiplies
2. Halving the RAM needed to store the coefficients in the SMP L2 cache

However, their technique is based on adjacent 4-choose-2 selectors (there being 6 possible combinations).

The problem with this is that there is an implicit assumption that a **dropout** that allows accurate results when
*random* items are adjacent (or semi-adjacent).  This may or may not be true. But it is a cheap trick used by
NVIDIA to inflate their performance 2..3x (2x being due to the reduction in actual multiplies (yet counting them in their marketing)
and another tiny fraction due to the decrease in L2 cache misses by halving the L2 cache memory pressure)).

A more robust solution would allow for `N-choose-N/2`, with a `dropout=0.5`.  Note, I doubt `0.5` is the best dropout.

Thus encoding a wider-scale bitmap that requires a full decoding step across a MUCH larger coefficient array should provide
equivalent reduction in multiplies along with a reduction in L2 memory pressure, while providing:
1. More flexible drop-out rates (need not adhere to exactly 50% dropout)
2. More flexible dropout positions (eg need not be semi-adjacent)

Further, I wish to explore the idea of using bit-set intersections instead of multiplications for **Attention** Networks.
I hypothesize that efficient 16-way attention selection (using 15 bitmaps with an average aggregate size of 4 bits per coefficient)
will achieve GREATER accuracy than current `qLORA` techniques (which quantize destructively to 4 bit coefficients).


# COMPARATIVE ANALYSIS

## FAX3
Early research utilized the **FAX3** Algorithm. This is an edge-detection algorithm used to transmit fonts over phone lines.
It works very well for having ranges of 1s, followed by ranges of 0s, and occasionally having blotches of 010101.
It is both fast and compresses VERY well.  However, it can not compress white-noise, and tends to grow such files.

## FAX4
A line-oriented previous-row comparison version of FAX3. If the preceeding line is a good predictor of the next line 
(which is true in FONTs in FAX), then each edge only needs about 2 bits.

## EWAH-32 and EWAH-64
Extended Word Aligned Hybrid.

This is a VERY high performance bitmap processing tool. In EWAH-32 You have a 15bit run-length followed by a 16bit literal length
and one bit to say if run is of 0s or 1s.  This is used in GIT for object-selection. You have a million objects, and the diff
says you only needs to grab 6 of them.. Then EWAH can describe that delta in a few bytes. By being word-aligned, you are not bit processing
and thus significantly faster than FAX.  EWAH works very well AND-combining and OR-combining pairs of bitmaps.

## 7/8 bit encoding of run-lengths with FAX4 line-pair diffs
The non word aligned flavor of the above EWAH. And the run-length version of FAX4.

Using an 8-bit number that uses the upper-8th bit as a need-more-bytes flag.  Thus 7bits is 1 byte, 14 bits is 2 bytes, etc.
This in combination with the FAX4 parent-row-diff technique means each edge within a 128-pixel boundary consumes exactly 1 byte.

If instead the gap distance between two line-edges is more than 128, it'll require 2 bytes per edge.

This, more-or-less is equivalent in speed/space as FAX4 and EWAH. Differing by being better in some situations, worse in others.

## Huffman
Many "fast" video processing frameworks just store the raw YCrCb (or RGB or RGBA) in a per-frame huffman table.
This allows FAX3 type speeds when decoding.

Huffman has the issue of LARGE table sizes. Thus for a 1024x1024 "frame", the overhead of the table is small, and thus worth it.
But for a small (how 512-choose-17) the table size would dwarf simply storing a 16bit array

## Bzip2, LZ77, LZW, LZ4, Snappy
These represent "dedup" encoding. For a bitmap (like a fax), it's unlikely to find repeated patterns, and thus these
techniques are not likely to be useful.  FAX3 or EWAH would do a better job of runs-of-zeros, which is the only realistic pattern.

## PDF jBIG2
PDF employs a range-encoding technique for a binary state of it's rasters. It uses a complex context to define a prediction for whether the next
pixel is white or black. That prediction produces a floating point probability.  It then uses a range-encoder to portray 
whether the prediction was correct or not.

As an example, if a lower-res image says this 8x8 is all black. The previous two rows all show all-black, and the pixel to your left (that has been decoded)
is black. Then you could say you have 99% probability of being black.. Lets codify that as a 2027 in 2048 chance of being black.

Thus we have a 11bit register and we ask, is it's value above 20. if so, then we emit black then SUBTRACT 20 from our integer.
If it crazily happened to be white (a 1 in 100 chance), then we have a value BELOW 20.. So we need to read 11 more bits in.

It's more complicated than this, and there is a dynamic component to the probabilities. But the essense is comparing a 11bit
number to a table of probabilities, comparing, then subtracting.

In my analysis, this performs very well for detected patterns - but the issue is in the just-barely-missed patterns.  Lets say
we are in a mode of predicting 001 001 001 type patterns. If a stretch instead encoded 0001100011, then we would be emitting 
MANY sets of 11 bits (to encode these 3 bit sequence fragments).

Similarly if I had an array of 1024 bits, in which exactly 1 was set. This would mis-predict with at least 11 bits across the whole
array.. In all likelihood this would require at least 23 bits to encode (range encoding typically wastest a full word worth of bits). 

Conversely an arithmatic approach would have to encode the length of N, the length of K.. Then it would be able to encode 1024-choose-1 with 10 bits.

My research therefore focused on situations where N-choose-K PLUS the size of the length-N and length-K would be similar-to 
or smaller-than the various flavors of range-encoding; as found in PDF jBIG2 and JPEG2000.  Even if the arithmatic flavor
doesn't compress as well, we can note that this 1024-choose-1 takes EXACTLY 3 computations.  length-N, length-K and read-log(N) 
v.s. 1024 macro calculation in jBIG (each requiring multiple steps)

## ANS (zstd)
On the heals of range-encoding comes AssymetricNumericalSystem.  This has two flavors: numbers and bits. The zstd which
is quickly replacing gzip and lzma for EVERYTHING is number-based. The ANS is the entropy portion of the deduping process.
It encodes huffman like sequences, but is more efficient when the distribution is not a perfect power-of-2 sequence.

Namely, shannon shows that you SHOULD be able to compress a probability distribution down to at least a certain level.
Huffman achieves this in the special case of perfect-powers-of-2 (example histograms: a=256,b=128,c=64). But this isn't
a practical histogram distribution.  ANS on the other hand nicely handles sub-ratios.. Typically used as a subdivision of 2048
sub-ranges (similar to range-encoding).

One interesting quirk; you have to decode in the opposite order as you encode (last byte to first), so it isn't as streamable
as gzip or lzma.

Another quirk is that, like range-encoding, you need a full 11 bit register. So your smallest encoding sizes are roughly 32bits
or so.  In practice, when coupled with a statistical model, you wind up with hundreds of bytes of overhead for the higher
fidelity histogram.

## Bernoulli Split - comparison

The reasons for continuing with BSplit given the playing field are partially due to engineering needs to produce compact 
regular-pattern or ultra-sparse bitmaps (1 random bit in 64).

At 1-delta-in-64 with NO pattern, range-encoding fails miserably (consuming at least 11 bits in each sub-sequence since it can't predict WHEN the bit will flip).

At 1-delta-in-64, pure ANS **can** do well. You store a single statistic (1/64), and it will semi-efficiently process all the bits.  
Howevever, you have to serially decode every bit.. So you'd need to perform 64 decodes to find the 1 bit.

BSplit, on the other hand ONLY encodes the minority bit.  So it runs the same if we have 1-in-64 or 63-in-64. Thus at these
ratios, our estimator above suggests `7.5 bits per element`. Just shy of a byte.  Note, this would be similar to
the 7/8 encoding of FAX4 above. HOWEVER, the 7/8 encoding can only get worse; e.g. if gaps are greater than 128, it uses 16 bits
per element.  BSplit, on the other hand is indifferent to WHERE the bits are; only that there is an average ratio of 1-to-64.
Thus we could have 200 1s bunched together with a gap of 12,800 and it would be (more or less) the same encoding size and
computational load.

## Heuristics

My research suggests the following average bit-size per sparseness ratio (ROUGH estimates)
```
avg-ratio  bits-per-element
1/2        1.999  (utilizing 41bit tail compaction)
1/3        2.999  (utilizing 41bit tail compaction)
1/4        3.75
1/5        4.0
1/8        4.5
1/10       5.0
1/16       5.5
1/24       6.0
1/64       7.5
1/96       8.0
1/128      8.5
1/192      9.0
1/256      9.5
```

Thus we can accurately predict the upper-bound compression of a bitset purely by determining the 1-to-0 ratio. Or
more specifically the minority / total.

Depending on the sparseness, the computation can skip over all the residuals when performing random access.

## Example

Example size 1000 items in a 64k address space.

```python
k=1000
n=64000
p=k/n                         # or 1 in 64
br=clg(n)-clg(k)              # 16-10 => 6 bits per element
resid_bits=bit_per_element*p  # 6*1000 => 6000 bits
pallet=2**clg(k+1)            # 1024 bits
collision_bits=k//2           # estimation of typical worst case; 512 bits (NOT a true upper bound)
bits=resid_bits+pallet+collision_bits # 6000+1024+512 => 7536 bits
avg_bits=bits/k               # 7536/1000 => 7.536 bits per element
predicted_bits_per_elem=7.5   # from table above or clg(64000)-clg(1000+1)+1.5 => 16-10+1.5 = 7.5
skip_bits(resid_bits,pallet,collision_bits,99) # skip 99 items, cursor will be on 100th item
```

The collision_bits could be slightly larger, but will typically be about 1/3rd of the pallet size for a random distribution.
Again, we are trading off performance for maximum compression. We could be more aggressive in extracting collision bits,
but that would remove some of our random access.

Thus to find the 100th item in a set of 1000 items (with an address space of 64 thousand possible items), we split
across many blocks of 64.. Each macro-block has between 1 and 64 items. We would read the pallet of each macro-block this
tells us the size of each collision-block.. Together we know how many bits to skip over (without actually decoding all the
individual elements).
Once in the target macro block (e.g. `k_left <= tgt_idx < k_left+ki`  ), then we skip over partition element-groups. We
need to decode the unitary-bit-vector to know how many elements are collided together.  If a majority of elements have no collisions,
then we can quickly skip over them without memory access.  We would need 1 memory access per collision skip detection.
For small br sizes (say br=6), most of those memory accesses would be in the same cache-line.  64*6 is only 48 bytes; thus
only needed 1 L1 cache miss.  An optimized block mapper could be used which assumes the entire residuals fits in a single
64bit register; thereby avoiding almost all per-item memory accesses.

This would be a preferred `CUDA` implementation (since memory access is the bottleneck).

# EXAMPLE USE CASES

* An alpha channel with a solid interior. Using a line-over-line predictor (such as FAX4) but with BSPlit, encoding densities can be nearly optimal AND decently fast.
* An array of 1-to-many sets (such as found in RDBMS's). Instead of encoding PrimaryKey-ForeignKey and ForeignKey-Primary key tuples. Store a compacted bitmap. 
This allows an in-memory database to be many orders of magnitude smaller than the traditional RDBMS model.
When coupled with a DataFrame library like Rust `Polars`. This can be a massive in-memory storage saving.
* A File system of offset+length pairs.  When storing inode to disk block mappings.. By storing all the block IDs in increasing order, 
and by trying to allocate blocks near each other (numerically), you can compact 32bit pointers into 8bit pointers (or even smaller).
* A classification search.  By categorizing entities into a tree and storing 1,2 or 3 BSplit bitmaps, you can "search" and scan
for items of class A or B or C.. You can also perform intersections and unions. If the bitmap is small enough to fit in L3 (or L2 or L1) cache
while the non-bitmap flavor would not; you can achieve higher search throughput. Though this is only possible if you tune the 
per-host entity-counts.
 