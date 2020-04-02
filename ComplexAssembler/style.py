
def get_writter(out_file):
	"""Returns a function that prints whatever string giben to it into
	the out_file can be used instead of print for easily redirect the messages
	to different files for different situations"""
	def fun(string):
		with open(out_file, "a") as myfile:
			myfile.write(str(string)+"\n")
	return fun
			
		
