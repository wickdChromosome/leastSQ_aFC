import logging

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
        """
        Call in a loop to create terminal progress bar
        @params:
                iteration   - Required  : current iteration (Int)
                total       - Required  : total iterations (Int)
                prefix      - Optional  : prefix string (Str)
                suffix      - Optional  : suffix string (Str)
                decimals    - Optional  : positive number of decimals in percent complete (Int)
                length      - Optional  : character length of bar (Int)
                fill        - Optional  : bar fill character (Str)
        """
        percent = str(iteration)
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print('\r%s |%s| %s %s' % (prefix, bar, percent, suffix), end = '\r')
        #Print New Line on Complete
        if iteration == total:
                print()


def intro(args):

	print("      ___         ___                                    ___                   ")
	print("     /\  \       /\__\        ___                       /\  \         _____    ")
	print("    /::\  \     /:/ _/_      /\__\                     /::\  \       /::\  \   ")
	print("   /:/\:\__\   /:/ /\__\    /:/__/                    /:/\:\  \     /:/\:\  \  ")
	print("  /:/ /:/  /  /:/ /:/ _/_  /::\  \     ___     ___   /:/ /::\  \   /:/ /::\__\ ")
	print(" /:/_/:/  /  /:/_/:/ /\__\ \/\:\  \   /\  \   /\__\ /:/_/:/\:\__\ /:/_/:/\:|__|")
	print(" \:\/:/  /   \:\/:/ /:/  /    \:\  \  \:\  \ /:/  / \:\/:/  \/__/ \:\/:/ /:/  /")
	print("  \::/__/     \::/_/:/  /      \:\__\  \:\  /:/  /   \::/__/       \::/_/:/  / ")
	print("   \:\  \      \:\/:/  /       /:/  /   \:\/:/  /     \:\  \        \:\/:/  /  ")
	print("    \:\__\      \::/  /       /:/  /     \::/  /       \:\__\        \::/  /   ")
	print("     \/__/       \/__/        \/__/       \/__/         \/__/         \/__/    ")
	print("										      ")
	print("										      ")
	print("					       					      ")
	print("AFC CALCULATOR from Pejlab - powered by ⅽ|_|				      ")
	print("Pejman Mohammadi | Bence Kotis						      ")
	print("										      ")
	print("#ARGUMENTS:							              ")
	print("#EQTL FILE: " + str(args.eqtl))
	print("#VCF/HAPLOTYPE FILE:  " + str(args.vcf))
	print("#EXPRESSIONS FILE: " + str(args.expr))
	print("#OUT FILE: " + str(args.output))
	print("										      ")
	print("										      ")

	logging.info("#ARGUMENTS:")
	logging.info("#EQTL FILE: " + str(args.eqtl))
	logging.info("#VCF/HAPLOTYPE FILE:  " + str(args.vcf))
	logging.info("#EXPRESSIONS FILE: " + str(args.expr))
	logging.info("#OUT FILE: " + str(args.output))
	logging.info("Rolled intro")
