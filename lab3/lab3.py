
import string
import matplotlib.pyplot as plt
f = open("romeoandjuliet.txt") #read the text
text = f.read()
f.close()

d = {}
words = []
punc = string.punctuation
content = text.lower().split()  

for w in content: #word processing
    if w[-1] in punc: # if ends with punc, remove
        new_w = w[:len(w)-1] 
        words.append(new_w)
    elif w[0] in punc: # if begins with punc, remove
        new_w = w[1:]
        words.append(new_w)
    else: # otherwise, add
        words.append(w)


for w in words: #add and count words to dictionary
    if w in d:
        d[w] += 1
    else:
        d[w] = 1

sorted_d = sorted(d.items(), reverse=True, key=lambda x:x[1]) #sorted in decreasing word count
top_ten = sorted_d[:10]
bottom_ten = sorted_d[-10:]

x = [x[0] for x in sorted_d]
y = [y[1] for y in sorted_d]

d_rank = {}
r = 1
for n in y:
    if n not in d_rank:
        d_rank[n] = r
        r += 1

n_words = [i for i in range(len(x))]  
plt.figure(num=None, figsize=(1000, 10), dpi=50)
plt.xlim(xmax=len(n_words))
plt.xlabel("words")
plt.ylabel("frequencies")
plt.bar(n_words,y)
plt.show()


plt.bar(x[:10], y[:10])
plt.show()

rank = [log(v) for v in d_rank.values()] 
count = [log(k) for k in d_rank.keys()]
plt.plot(count, rank)
plt.xlabel('word count')
plt.ylabel('word rank')
plt.show()

