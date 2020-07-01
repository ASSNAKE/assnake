
So

check for config

manage all params.json or .yaml files
manage all static files


different name different content
same name and content, hash not written to filename
same content different names
same name different contents

a 1
b 2
b 1




Say, we have `good` params in dataset and `good` in database but the contents are different.
We are importing dataset in assnake - collision!
We defenitely need `hash`, but how can we elegantly manage it in workflow and resulting folder and filenames?

```
input:  params  = '{db}/{tool}/{params_name}.json'
output: primary = './raw__tool_{params_name}/{---}.out'
```
`dataset/assnake_data/params/{tool}/{param_name}.json`

Or let's say that one unintentionally modifies `params.json` and then tries to run the tool again. 
One will get an error message that this parameters are different. 

So, okay what if we carry CRC32 hash all the time:
```
input:  params  = '{db}/{tool}/{params_name}:{params_hash}.json'
output: primary = './raw__tool_{params_name}:{params_hash}/{---}.out'
```

It makes paths longer and more unreadable, and certainly untypable, but it only really matters in cases of basename collision. 
Otherwise we should make the user believe that there is no hashes at all. 

