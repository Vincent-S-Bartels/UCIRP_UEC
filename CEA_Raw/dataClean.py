def cleanData(inFile, outFile, header = 0):
    with open(inFile, 'r') as f1:
        x = f1.readlines()
        x = x[header:]
        with open(outFile, 'w') as f2:
            for i in x:
                f2.write(','.join(i.split()) + '\n')
