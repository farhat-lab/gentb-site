


def unpack_mutation_format(mutation):
    """
    Processes a mutation into it's three main parts:

     * Index - Optional order of the mutation
     * Gene - The gene involved in this mutation
     * Mutation - Returned Mutation, usually just cleaned.

    """
    index = None
    if " " in mutation:
	index, mutation = mutation.split(" ", 1)
	try:
	    index = int(index)
	except:
	    raise ValueError("Optional sort index should be a number.")

    bits = mutation.split('_')
    if bits[1] == 'P':
	if 'promoter' not in bits:
	    raise ValueError("Promoter doesn't specify 'promoter' part")
	name = ' '.join(bits[bits.index('promoter') + 1:])
	return (index, 'promoter ' + name, mutation)
    elif bits[1] == 'I':
	if 'inter' not in bits:
	    raise ValueError("Integenic doesn't specify 'inter' part")
	name = ' '.join(bits[bits.index('inter') + 1:])
	return (index, 'intergenic ' + name, mutation)
    elif bits[1] in ['CN', 'CD', 'CF', 'CI', 'CZ', 'N', 'ND', 'NI', 'NF']:
	return (index, bits[-1], mutation)

    raise ValueError("Must be promoter, intergenic or CN, CD, CF, CI, CZ or N, ND, NI, NF")

