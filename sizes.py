import unittest
import math


def lg(k: int) -> int:
    """log_2(k)"""
    return int(math.log2(k))

def clg(k: int) -> int:
    """ceil(log_2(k))"""
    return int(math.ceil(math.log2(k)))
def flg(k: int) -> int:
    """floor(log_2(k))"""
    return int(math.floor(math.log2(k)))

def cflg(k: int)->(int,int):
    """Calculate BOTH: ceil(log_2(k)), floor(log_2(k))"""
    x=math.log2(k)
    return (int(math.floor(x)),int(math.ceil(x)))

def nck(n: int, k: int) -> float:
    """log_2(N-Choose-K)"""
    return math.log2(math.comb(n, k))


# global cache - NOT for production use
cnk_cache = {}


def cnk(n: int, k: int) -> int:
    """Cached ceil(log_2(N-Choose-K))
    Be wary of using this in production for unbounded N,K; you will run out of RAM"""
    global cnk_cache
    val = cnk_cache.get((n, k), None)
    if val is None:
        val = int(math.ceil(nck(n, k)))
        cnk_cache[(n, k)] = val
    return val


def lr_nck_key_size(k, k_left) -> int:
    """Size of the finite state machine for the left/right split"""
    if k == 3:
        return 2
    if k == 2:
        if k_left in {0, 2}:
            return 2
        else:
            return 1
    if k == 5:
        if k_left in {2, 3}:
            return 2
        else:
            return 3
    if k == 6:
        if k_left in {2, 3, 4}:
            return 2
        else:
            return 4
    if k == 8:
        if k_left in {3, 4, 5}:
            return 2
        elif k_left in {2, 6}:
            return 4
        else:
            return 5
    if k == 10:
        if k_left in {4, 5, 6}:
            return 2
        elif k_left in {3, 7}:
            return 4
        elif k_left in {2, 8}:
            return 5
        else:
            return 6
    k2 = k // 2
    delta = k_left - k2 - 1 if k_left > k2 else k2 - k_left
    if k & 0b1 == 0:
        # even
        if k <= 6 or k == 10:
            if delta < 2:
                return 2
            elif k_left == 0 or k_left == k:
                return k2 + 1
            else:
                return 4 + delta - 2
        else:
            if delta < 4:
                return 3
            elif delta == k2:
                return k2
            else:
                return 1 + delta
    else:
        # odd
        if delta == k2:
            return 1 + delta
        else:
            return 2 + delta


# globals for use with internal investigation of function ranges
cache = {}
worsts = []
worsts_table = {}


def show_worsts():
    """Export worst-case sizes"""
    return [(k, worsts_table[k]) for k in sorted(worsts_table)]


def clear_worsts():
    """Reset global statistics"""
    global cache
    cache = {}
    global worsts
    worsts = []
    global worsts_table
    worsts_table = {}


# Tolerance hack - TODO remove
hack = False


def bm_size(n) -> int:
    """Size of the BitMap for a given N
    For n:1..8 we just use pairty to determine last bit (thus n-1)
    For n:9..40 we emit the leading bits then calculate an arithmatic 9-choose-K
    Currently using a 3*2**4 as it's more convenient, thus 9..47
    either way, 9-choose-4 is worst case; 7 bits.. 9-7 leaves 2 bits saved
    For else: emit prefix bits, then use arithmatic 41-choose-K (currently 48-choose-K)
    either way, 41-choose-20 is 38, 48-choose-24 is 45; both save 3 bits.
    The only difference is if N is between 41 and 48, the 41-choose-K would be 1 bit smaller
    """
    if n < 9: return n - 1
    if n < 48: return n - 2
    return n - 3


def lr_nck_max(n, k) -> int:
    """Try and find the max size (exhaustively) of N-Choose-K using the left/right split, and our cache"""
    global cache
    global worsts
    key = (n, k)
    if key in cache: return cache[key]
    n2 = n // 2
    n_right = n - n2
    if k == 0:
        val = 0
    elif k == 1:
        val = clg(n)
    elif k > n2:
        val = lr_nck_max(n, n - k)
    elif 3 * k > n:  # assumes arithmatic
        return bm_size(n)
    elif hack and n == 16 and k == 4:
        return 11  # HACK!!!!!!
    else:
        max1 = 0
        max_tup = None
        _cnk = cnk(n, k)
        for ki in range(0, k + 1):
            key_size = lr_nck_key_size(k, ki)
            probe = lr_nck_max(n2, ki) + lr_nck_max(n_right, k - ki) + key_size
            if probe > max1:
                max1 = probe
                max_tup = ((n, k, ki), key_size, probe, _cnk)
        val = max1
        if not max_tup is None:
            worsts.append(max_tup)
            worsts_table[f'{n:02}-{k:02}-{max_tup[0][2]:02}'] = max_tup[1:]
    val = int(val)
    cache[key] = val
    return val


def lr_nck_est(n, k) -> int:
    """Approximate the bit-size of N-choose-K using the left/right split"""
    if k <= 8:
        if k == 9:
            r = 16
        elif k == 8:
            r = 13
        elif k == 7:
            r = 11
        elif k == 6:
            r = 8
        elif k == 5:
            r = 6
        elif k == 4:
            r = 4
        elif k == 3:
            r = 2
        elif k == 2:
            r = 1
        return clg(n) * k - r
    return 0


def lr_nck_max_bn(bn, k) -> int:
    """Calculate worst-case size of left/right split"""
    n = 2 ** bn
    return lr_nck_max(n, k)


def inner_hat(bn: int, k: int, d: int, inner) -> int | None:
    """This is the more stable larger K-size split (16<=K<=64).
    We use a pallet of 2**bk to partition bn into approximately K partitions of size 2**bk each.
    We use a collision pallet (with multiple encoding methods) to determine if adjacent br entries are convolved.
    We use a unitary bit-sequence to identify the size of the collision..
    So "1" is length 1 (no collisions), br + 1 bit
       "01" is length 2 (1 collision - allows 1 bit to be extracted with lr_split; net size is br + 1 bits)
       "001" is length 3 (2 collisions - allows 2 bits to be extracted with lr_split; net size is br + 1 bits)
       "0001" is length 4 (3 collisions - allwos 4 bits to be extracted with lr_split; net size is br bits)
       and so on, larger collisions reduce the net size by larger amounts; however, to allow for random-access, we
       must SKIP-OVER those unmapped bits.  We thus prefer randomized keys that minimize the probability of collisions.
    By storing the first bit of the collision sequence for a given element in a separate map, we can maintain a precise
    br*k array. Since for the 1,2,3 collision-length cases, we have the exact same residual size.

    Note, if we have more than 8 collisions, we can recurse inner_hat, since we haven't defined lr_split for k > 8
    """
    n = 2 ** bn
    if k == 0:
        return 0
    if k > n:
        return None
    hn = n // 2
    if k > hn:
        return inner_hat(bn, n - k, (n - k) // 2, inner)
    if k * 3 > n:
        return bm_size(n)  # client/server know this
    if d >= k:
        raise ValueError(f'd must be less than k: {d} >= {k}')
    u = k - d
    if u == 0:
        raise ValueError(f'u must be greater than 0: {u}')
    bk = clg(k + 1)
    pallet = 2 ** bk
    if u == 1:
        sig = inner_hat(bn - bk, k, 0, inner)
        if sig is None:
            return None
        return sig + pallet
    if d == 0:
        advance = 0
    else:
        min_val = clg(k - 1) * (u - 1) - d
        if d > 0:
            chk_cbu = clg(u) * (d - 1) - d
            if chk_cbu < min_val:
                min_val = chk_cbu
        if min_val > u - 1:
            min_val = u - 1
        half_k = (k + 1) // 2 + 8
        if min_val <= half_k:
            advance = min_val
        else:
            advance = inner(k - 1, d - 1) - d
    resid = (bn - bk) * k
    tot = pallet + advance + resid
    if tot > n - 3:
        return n - 3
    return tot


def hat_a(bn: int, k: int, d: int) -> int | None:
    """Investigative - use arithmatic sizing to calculate inner bit-size"""
    return inner_hat(bn, k, d, lambda a, b: cnk(a, b))


def hat_lr(bn: int, k: int, d: int) -> int | None:
    """Investigative - use left/right split sizing to calculate inner bit-size"""
    return inner_hat(bn, k, d, lambda a, b: lr_nck_max(a, b))


def hat_est(bn, k) -> int:
    if bn < 2: raise ValueError(f'bn must be at least 2: {bn}')
    n = 2 ** bn
    hn = n // 2
    if k < 0: raise ValueError(f'k must be non-negative: {k}')
    if k == 0: return 0
    if k == 1: return clg(n)
    if k > hn:
        return hat_est(bn, n - k)
    if k * 3 > n:
        return bm_size(n)
    bk = clg(k)
    pallet = 2 ** bk
    advance = (k + 1) // 2
    resid = (bn - bk) * k
    r = pallet + advance + resid
    alt_size = bm_size(n)
    if r > alt_size:
        return alt_size + 1
    return r


a_41_v = [nck(38, i) for i in range(0, 39)]
a_38_v = [nck(35, i) for i in range(0, 36)]
a_35_v = [nck(32, i) for i in range(0, 33)]
a_32_v = [nck(29, i) for i in range(0, 30)]
a_29_v = [nck(26, i) for i in range(0, 27)]
a_26_v = [nck(23, i) for i in range(0, 24)]
a_23_v = [nck(20, i) for i in range(0, 21)]
a_20_v = [nck(17, i) for i in range(0, 18)]
a_17_v = [nck(14, i) for i in range(0, 15)]
a_14_v = [nck(11, i) for i in range(0, 12)]
a_11_v = [nck(9, i) for i in range(0, 10)]
a_9_v = [nck(6, i) for i in range(0, 7)]
a_6_v = [nck(3, i) for i in range(0, 4)]

a_t = [0 for _ in range(0, 42)]
a_t[6] = a_6_v
a_t[9] = a_9_v
a_t[11] = a_11_v
a_t[14] = a_14_v
a_t[17] = a_17_v
a_t[20] = a_20_v
a_t[23] = a_23_v
a_t[26] = a_26_v
a_t[29] = a_29_v
a_t[32] = a_32_v
a_t[35] = a_35_v
a_t[38] = a_38_v
a_t[41] = a_41_v
a_n = [3 for _ in range(0, 42)]
a_n[1] = 2


def trial_a_x(n, k: int, bits: int) -> int:
    if k == 0: return 0
    if 2 * k > n:
        mask = 2 ** n - 1
        return (~trial_a_x(n, n - k, bits)) & mask
    if k == 1: return 1 << bits
    num_left_bits = a_n[n]
    mask = 2 ** num_left_bits - 1
    left_bits = bits & mask
    k_left = left_bits.bit_count()
    k_right = k - k_left
    bits = bits >> num_left_bits
    n_nxt = n - num_left_bits
    table_this = a_t[n]
    table_nxt = a_t[n_nxt]
    if num_left_bits == 3:
        if k_left == 0:
            # 1 permutation
            nxt_k = k
            nxt_bits = bits
        elif k_left == 1:
            # 3 permutations
            nxt_k = k - 1
            nxt_bits = bits - table_nxt[k_right]
        elif k_left == 2:
            # 3 permutations
            nxt_k = k - 2
            nxt_bits = bits - table_nxt[k_right] - table_nxt[k_right - 1] * 3
        else:
            # k_left == 3 1 permutation
            nxt_k = k - 3
            nxt_bits = bits - (table_this[k] - table_nxt[k - 3])
    else:
        if k_left == 0:
            # 1 permutation
            nxt_k = k
            nxt_bits = bits
        elif k_left == 1:
            # 2 permutations
            nxt_k = k - 1
            nxt_bits = bits - table_nxt[k_right]
        else:
            # k_left == 2  1 permutation
            nxt_k = k - 2
            nxt_bits = bits - (table_this[k] - table_nxt[k - 2])
    return left_bits | (trial_a_x(n_nxt, nxt_k, nxt_bits) << num_left_bits)


def trial_a_41_plus(n: int, k: int, bits: int) -> int:
    left = n - 41
    mask = 2 ** left - 1
    out = bits & mask
    kl = out.bit_count()
    kr = k - kl
    bits = bits >> left
    return out | (trial_a_x(kr, bits) << left)


def trial_a_8(k: int, bits: int) -> int:
    left = n - 9
    mask = 2 ** left - 1
    out = bits & mask
    kl = out.bit_count()
    kr = k - kl
    bits = bits >> left
    return out | (trial_a_x(kr, bits) << left)


def trial_a_9_plus(n: int, k: int, bits: int) -> int:
    left = n - 9
    mask = 2 ** left - 1
    out = bits & mask
    kl = out.bit_count()
    kr = k - kl
    bits = bits >> left
    return out | (trial_a_x(kr, bits) << left)


def trial_a(n: int, k: int, bits: int) -> int:
    if n < 41:
        if n < 9:
            if n == 8: return trial_a_8(k, bits)
            if n == 7: return trial_a_7(k, bits)
            if n == 6: return trial_a_6(k, bits)
            if n == 5: return trial_a_5(k, bits)
            if n == 4: return trial_a_4(k, bits)
            if n == 3: return trial_a_3(k, bits)
            if n == 2: return trial_a_2(k, bits)
            return trial_a_1(k, bits)
        else:
            return trial_a_9_plus(n, k, bits)
    else:
        return trial_a_41_plus(n, k, bits)


def reportse(bn, s, e):
    n = 2 ** bn
    for k in range(s, e):
        print()
        h_e = hat_est(bn, k)
        lr = lr_nck_max_bn(bn, k)
        ccnk = cnk(n, k)
        arr = clg(n) * k
        min1 = 9999999
        max1 = 0
        first = 0
        for d in range(0, k):
            h_a = hat_a(bn, k, d)
            if h_a is None:
                continue
            h_lr = hat_lr(bn, k, d)
            if h_lr is None:
                continue
            if d == 0:
                first = h_lr
            if h_lr > h_e and h_a > h_e:
                tok = 'LR A'
            elif h_lr > h_e:
                tok = 'LR  '
            elif h_a > h_e:
                tok = '   A'
            else:
                tok = '    '
            print(
                f'{tok} k:{k:03} d:{d:03} h_lr:{h_lr:04} h_a:{h_a:04} (he:{h_e:04} lr:{lr:04} nk:{ccnk:04} arr:{arr:04})')
            min1 = min(min1, h_lr)
            max1 = max(max1, h_lr)
        spread = max1 - min1
        spread1 = max1 - first
        print(f'---- k:{k:02} min:{min1:03} max:{max1:03} spread:{spread:03} spread0:{spread1:03}')


def report1(bn, k):
    reportse(bn, k, k + 1)


def q_size(v:int, max:int)->int:
    bnm1,bnp1 = cflg(max)
    T=2**bnp1-max
    if v < T:
        return bnm1
    return bnp1

def report():
    reportse(8, 10, 40)

def jjj_suffix(k: int, u: int, d: int, t: int = 33) -> int:
    if d == 0:
        bits = 0
    elif d * clg(u) < u:
        bits = d * clg(u)
    elif u < d:
        bits = u
    else:
        bits = min(u,clg(u) + cnk(u, d))
    return bits


def jjj(min1: int, max1: int = None, show: bool = True):
    if max1 == None:
        max1 = min1 + 1
    worst = 0
    for k in range(min1, max1):
        bias = 0
        if k < 32:
            if k < 16:
                bias = k * 2
                v = 16
            else:
                bias = k
                v = 32
        else:
            v = 64
        t = int(v / 2) + 1
        clv = clg(v)
        for u in range(1, k + 1):
            d = int(k - u)
            bits = clv + cnk(v, u) + bias
            bits += jjj_suffix(k, u, d, t)
            if show:
                print(f'k:{k:02} u:{u:02} d:{d:02} bits:{bits:03}')
            worst = max(worst, bits)
    if show:
        print(f'worst:{worst:03}')
    return worst


def lll(min1: int, max1: int = None, show: bool = True):
    if max1 == None:
        max1 = min1 + 1
    worst = 0
    for k in range(min1, max1):
        bias = 0
        if k < 32:
            if k < 16:
                bias = k * 2
                v = 16
            else:
                bias = k
                v = 32
        else:
            v = 64
        t = int(v / 2) + 1
        for u in range(1, k + 1):
            d = int(k - u)
            bits = v + bias
            bits += jjj_suffix(k, u, d, t)
            if show:
                print(f'k:{k:02} u:{u:02} d:{d:02} bits:{bits:03}')
            worst = max(worst, bits)
    if show:
        print(f'worst:{worst:03}')
    return worst


def foo():
    for k in range(10, 65):
        bits=cnk(2**16,k)-(16-6)*k
        print(f'k:{k:03}: {jjj(k, show=False)} {lll(k, show=False)} {bits}')

def trial_sum(name:str,n:int,k:int)->str:
    sum=0
    name = f'{name}_{n}_{k}';
    half_n = n//2
    half_k = k//2
    end = min(k, n // 2)
    total=0
    portion=0
    bn=clg(n)
    if k%2==0:
        # even
        q = math.comb(half_n,half_k)
        v = q**2
        sum=v
        total=v
        txt = f'const {name}: WE<{end},{bn}> = WE {{ even: CdfQ {{ cdf: {v}, q: {q} }}, table: ['
    else:
        # odd
        txt = f'const {name}: WO<{end},{bn}> = WO {{ table: ['
    for i in range(0,end):
        sz=math.comb(half_n,i)
        tsz=sz * math.comb(half_n,k-i)
        sum += tsz
        portion += sz
        txt = f'{txt} CdfQ {{ cdf: {sum}, q: {sz} }},'
    total += portion * 2
    if total != math.comb(half_n,k):
        print("ERROR")
    txt += f' ], total: {total} }};'
    return txt

import random


def binary_search(arr, x)->(bool,int):
    low = 0
    high = len(arr) - 1
    mid = 0
    while low <= high:
        mid = (high + low) // 2
        if arr[mid] < x:
            low = mid + 1
        elif arr[mid] > x:
            high = mid - 1
        else:
            return (True,mid)
    return (False,low)

def gen_rand(n:int)->list[int]:
    return [random.randint(0,2**32-1) for _ in range(0,n)]

def gen_rand_max(n:int,max:int)->list[int]:
    mm1=max-1
    return [random.randint(0,mm1) for _ in range(0,n)]

def shrink_right(items:list[int], n_left:int)->list[int]:
    return [x-n_left for x in items]

def invert(items:list[int],max:int)->list[int]:
    if len(items) == 0:
        return [x for x in range(0,max)]
    if len(items) == 1:
        return [x for x in range(0,items[0])] + [x for x in range(items[0]+1,max)]
    if len(items) == 2:
        return [x for x in range(0, items[0])] + [x for x in range(items[0] + 1, items[1])] + [x for x in range(items[1] + 1, max)]
    res=[]
    next_v=items[0]
    last_i=0
    items_i=1
    while True:
        for x in range(last_i, next_v):
            res.append(x)
        last_i = next_v + 1
        if items_i >= len(items):
            break

        next_v=items[items_i]
        items_i += 1
    for x in range(last_i, max):
        res.append(x)
    return res

def two_depth(items:list[int],bn:int)->int:
    if bn <= 1:
        pass
    assert bn > 1
    depth=0
    k=len(items)
    while True:
        n=2**(bn-1)
        if n <= k:
            return depth
        if items[0] < n:
            if items[-1] >= n:
                return depth
        else:
            return depth
        depth += 1
        bn -= 1
        if bn <= 1:
            pass
        assert bn > 1

def route_size(items:list[int],bn:int,max_bn:int=None)->int:
    if max_bn is None:
        max_bn=bn
    k=len(items)
    if k == 1 or k == 0:
        return 0
    assert k > 1
    assert bn > 1
    n=2**bn
    if 2*k > n:
        print("invert")
        return route_size(invert(items,2**bn),bn, max_bn)
    if 2*k == n:
        print('mid')
        return n
    # if 3*k > n:
    #     print(f"Too Big n:{n} k:{k} bn:{bn}")
    #     return n
    bk_floor=flg(k)
    td=two_depth(items,bn-bk_floor)
    bits=q_size(td,bn-bk_floor)
    half=2**(bn-1)
    f,p=binary_search(items,half)
    left=items[0:p]
    right=shrink_right(items[p:],half)
    k_left=len(left)
    k_right=k-k_left
    assert k_right == len(right)
    pad=' '*(max_bn-bn)
    # print(f'{pad}[{k_left}/{k_right}]')
    bn_nxt=bn-1-td
    bits += route_size(left,bn_nxt,max_bn)
    bits += route_size(right,bn_nxt,max_bn)
    # print(f'{pad}<{k_left}/{k_right}>: {bits}')
    return bits

class Test(unittest.TestCase):
    def test_search(self):
        nums=gen_rand(100)
        nums.sort()
        assert (True,0) == binary_search(nums,nums[0])
        assert (True, len(nums)-1) == binary_search(nums, nums[-1])
        assert (False,1) == binary_search(nums,nums[1]-1)

    def test_1(self):
        foo()
    def test_2(self):
        v=cnk(5,5)

    def test_3(self):
        print(f'{invert([2], 10)}')
        assert [0, 1, 3, 4, 5, 6, 7, 8, 9]==invert([2], 10)
        print(f'{invert([2,4], 10)}')
        assert [0, 1, 3, 5, 6, 7, 8, 9]==invert([2,4], 10)
        print(f'{invert([2, 4, 6], 10)}')
        assert [0, 1, 3, 5, 7, 8, 9] == invert([2, 4, 6], 10)
        print(f'{invert([2,4,6,8],10)}')
        assert [0, 1, 3, 5, 7, 9] == invert([2,4,6,8],10)
        print(f'{invert([0,1,2,3], 4)}')
        assert [] == invert([0,1,2,3], 4)

    def test_old_main(self):
        report()
        hat_a(24, 6, 2)

    def test_q(self):
        for j in [2,3,4,6,8,12,14]:
            for i in range(0, j):
                print(f'q_size({i},{j}): {q_size(i,j)}')


    def test_search(self):
        bn=9
        nums=gen_rand_max(255,2**bn)
        nums.sort()
        print(f'sz: {route_size(nums,bn)}')
