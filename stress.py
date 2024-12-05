import sys

def consume_memory(n):
    memory = []
    while True:
        memory.append(' ' * 10**6)  # Append a string that takes ~1 MB
        print(f"Allocated {len(memory)} MB")
        if len(memory) >= n:
            break
          
mem = sys.argv[1]  # Consume ~1 GB of memory
consume_memory(int(mem) * 1024)  
