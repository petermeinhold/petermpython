[3:40:16 AM] Andrea Zonca: sure
[3:42:10 AM] Andrea Zonca: found linear bin:
[3:42:12 AM] Andrea Zonca: def bin(y, samples):
    y_split = np.array_split(y, len(y)/samples)
    return map(np.mean, y_split), np.array(map(np.std, y_split))/np.sqrt(samples)
[3:42:18 AM] Andrea Zonca: looking for logbin now
[3:43:03 AM] Peter Meinhold: hmm
[3:43:08 AM] Peter Meinhold: is this in python?
[3:43:27 AM] Andrea Zonca: yes
[3:43:34 AM] Andrea Zonca: let's say samples is 5
[3:43:50 AM] Andrea Zonca: I split the array in 5 samples chunks
[3:43:53 AM] Andrea Zonca: with array_split
[3:44:04 AM] Andrea Zonca: then I apply to each chunck the np.mean operator
[3:44:12 AM] Andrea Zonca: map ( function, list)
[3:44:21 AM] Andrea Zonca: applies the function to each element of a list