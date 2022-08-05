import functools
from operator import attrgetter


class Person:
    def __init__(self,name,work,score,loss):
        self.name = name
        self.work = work
        self.score = score
        self.loss = loss

    def __repr__(self):
        return self.name




def compare(x,y):
    if x.score == y.score:
        if x.loss < y.loss:
            return True
        return False
    else:
        return x.score > y.score


def partition(arr, l, r):
    x = arr[r]
    i = l
    for j in range(l, r):
        if compare(arr[j], x):
            arr[i], arr[j] = arr[j], arr[i]
            i += 1
    arr[i], arr[r] = arr[r], arr[i]
    return i


def kth(arr, l, r, k):
    if (k > r - l + 1):
        return 
    index = partition(arr, l, r)
    if (index - l == k - 1):
        return arr[index]
    if (index - l > k - 1):
        return kth(arr, l, index - 1, k)
    return kth(arr, index + 1, r,
                       k - index + l - 1)


    

n = int(input())

accepted_candidates = []
index_name = [0] * n
works = { }
max_workers = { }
for i in range (n):
    s, m = input().split(",")
    m = int(m)
    works[s] = []
    max_workers[s] = m

k = int(input())

for i in range(k):
    c,q,r,p = input().split(",")
    r = int(r)
    p = int(p)
    person = Person(c,q,r,p)
    works[q].append(person)

for key, value in works.items():
    max = max_workers[key]
    if max < len(max_workers[key]):
        kth(works[key], 0, len(max_workers[key]), max)
    print(works[key])
    for i in range(min(max,len(works[key]))):
        accepted_candidates.append(works[key][i].name)

for x in sorted(accepted_candidates):
    print(x)